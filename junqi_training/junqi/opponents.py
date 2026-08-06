"""The opponent pool a candidate is judged against.

Beating only the immediately previous version is how a bot cycles: it learns to
exploit one opponent and forgets everything else. Every candidate here plays the
whole pool, and its fitness is the aggregate.

Agents are described by a picklable :class:`AgentSpec` and built inside the
worker process, so a pool can be fanned out across cores.
"""

from __future__ import annotations

import random
from dataclasses import dataclass, field, replace
from importlib import import_module
from pathlib import Path

from .bot import (
    PRIOR_BATTLE,
    BotWeights,
    HeuristicBot,
    _distance,
    _expected_battle,
    _piece_value,
    enemy_flag_squares,
)
from .deployment import rear_rows
from .game import Game
from .knowledge import OpponentKnowledge
from .search_bot import SearchBot
from .types import Move, Owner, PieceKind, Position


class RandomBot:
    """Floor of the pool. A candidate that cannot beat this is broken."""

    def __init__(self, seed: int | None = None) -> None:
        self.rng = random.Random(seed)

    def choose_move(self, game: Game, owner: Owner | None = None) -> Move:
        moves = game.legal_moves(game.turn if owner is None else owner)
        if not moves:
            raise ValueError("当前玩家没有合法走法")
        return self.rng.choice(moves)


class HQRushBot:
    """The strategy a human finds first: walk cheap pieces at the enemy flag.

    This is the opponent the bot was actually losing to, so it belongs in the
    pool permanently rather than as a one-off benchmark.
    """

    def __init__(self, seed: int | None = None, caution: float = 1.0) -> None:
        self.rng = random.Random(seed)
        self.caution = caution

    def choose_move(self, game: Game, owner: Owner | None = None) -> Move:
        player = game.turn if owner is None else owner
        targets = game.flag_candidates(player.other)
        moves = game.legal_moves(player)
        if not moves:
            raise ValueError("当前玩家没有合法走法")
        best: Move | None = None
        best_score = None
        rear = {0, 1} if player == Owner.HUMAN else {10, 11}
        for move in moves:
            piece = game.board[move.src]
            score = self.rng.uniform(-0.4, 0.4)
            if targets:
                before = min(_distance(move.src, square) for square in targets)
                after = min(_distance(move.dst, square) for square in targets)
                score += (before - after) * 10.0
                if move.dst in targets:
                    score += 1_000.0
            score -= _piece_value(piece.kind) * 0.5
            occupant = game.board.get(move.dst)
            if occupant is not None and move.dst not in targets:
                score -= 4.0 * self.caution
                if move.dst[0] in rear and piece.kind != PieceKind.ENGINEER:
                    score -= 8.0 * self.caution
            if best_score is None or score > best_score:
                best_score, best = score, move
        assert best is not None
        return best


class OracleSearchBot(SearchBot):
    """A deliberately cheating opponent: it reads every hidden rank.

    **Never ships, never plays a human.** It exists because the judging pool has
    no ceiling: `search-mid` tops out at about 60% against the subject and *is*
    the subject's own policy class carrying anchor weights, so the pool cannot
    grade anything past parity-with-self. An omniscient opponent is structurally
    different in the one dimension this game is about, and it is free to build.

    This does not violate the hidden-rank invariant. That invariant binds the
    *subject* -- the agent that ships. Using privileged information in a
    training or evaluation opponent is standard practice; Suphx's oracle agent
    is the same device.

    Measured at equal weights and equal search widths, the oracle side wins
    0.604 +/- 0.023 over 400 games, so perfect information is worth about +10
    points of win rate here. That makes it hard but beatable -- the right shape
    for a pool ceiling rather than an unwinnable wall.

    `gamma` below 1 reveals only that fraction of the hidden pieces, which is
    Suphx's dropout schedule and gives a dial between the belief bot and full
    clairvoyance.
    """

    def __init__(
        self,
        weights: BotWeights | None = None,
        seed: int | None = None,
        samples: int = 6,
        beam_width: int = 10,
        reply_width: int = 4,
        gamma: float = 1.0,
        gamma_beam: float = 0.0,
    ) -> None:
        # At gamma 1 every sampled world is the same world, so N samples is N
        # identical rollouts. Collapsing to one makes the strongest opponent in
        # the pool also the cheapest.
        super().__init__(
            weights,
            seed,
            1 if gamma >= 1.0 else samples,
            beam_width,
            reply_width,
        )
        self.gamma = gamma
        self.gamma_beam = gamma_beam

    def _update_knowledge(self, game: Game, owner: Owner) -> None:
        """Optionally let the *candidate generator* cheat as well.

        `gamma` alone only buys a truthful rollout. The beam that decides which
        ten moves get rolled out at all is still ranked by
        `HeuristicBot._score` using ordinary deduced belief -- so a move that is
        strong *only because* you can see the enemy's ranks ("take C7, it is a
        lieutenant") has to look attractive to a blind heuristic before the
        oracle is even allowed to consider it. That caps how much the cheat is
        worth, and it is why the plain oracle understates perfect information.

        Setting `gamma_beam` collapses the belief to singletons of the true
        ranks, so the heuristic prices every capture exactly and
        `enemy_flag_squares` resolves to the real flag square. This is the
        genuine upper bound, and correspondingly the least like a human.
        """
        super()._update_knowledge(game, owner)
        if self.gamma_beam <= 0.0 or self.knowledge is None:
            return
        # Seeded from the ply so a given position always reveals the same
        # subset, keeping the agent deterministic under common random numbers.
        rng = random.Random(self.base_seed * 7_919 + game.move_count)
        for position, piece in game.board.items():
            if piece.owner != owner.other or piece.revealed:
                continue
            if self.gamma_beam >= 1.0 or rng.random() < self.gamma_beam:
                self.knowledge.possible[position] = frozenset({piece.kind})

    def _determinize(
        self, game: Game, hidden_owner: Owner, rng: random.Random
    ) -> Game:
        if self.gamma >= 1.0:
            return game.clone()
        sampled = super()._determinize(game, hidden_owner, rng)
        for position, piece in game.board.items():
            if (
                piece.owner == hidden_owner
                and not piece.revealed
                and rng.random() < self.gamma
            ):
                sampled.board[position] = piece
        return sampled


class SelectiveBot:
    """Picks its fights: attacks only when the odds are with it.

    The policy class the pool was missing, and the reason a mispriced attack
    term could survive training. Every weight-driven opponent here shares
    `capture` and `blind_battle` with the candidate -- `material_weights` even
    doubles `capture` -- so the whole pool over-attacks in the same direction
    and the bias nets to zero in self-play. `HQRushBot` is the only
    structurally different opponent and it carries the *opposite* bias, since
    it refuses captures outright. Nothing attacked *selectively*, which is what
    a human does, and over ten replayed games the human won 76% of the battles
    it started against the bot's 46%.

    Deduces from its own battle history exactly as `SearchBot` does, so it
    never reads a hidden rank -- it just declines fights a prior says it loses.
    """

    def __init__(
        self,
        weights: BotWeights,
        seed: int | None = None,
        threshold: float = 0.05,
    ) -> None:
        self.heuristic = HeuristicBot(weights, seed=seed)
        self.rng = random.Random(seed)
        self.threshold = threshold
        self.knowledge: OpponentKnowledge | None = None
        self.processed_records = 0

    def choose_move(self, game: Game, owner: Owner | None = None) -> Move:
        player = game.turn if owner is None else owner
        self._update_knowledge(game, player)
        assert self.knowledge is not None
        self.heuristic.knowledge = self.knowledge.possible
        moves = game.legal_moves(player)
        if not moves:
            raise ValueError("当前玩家没有合法走法")
        targets = enemy_flag_squares(game, player, self.knowledge.possible)
        # Fall back to the full list when every move is a fight it dislikes --
        # declining to move is not an option.
        allowed = [
            move for move in moves if self._acceptable(game, move, player, targets)
        ] or moves
        scored = [
            (self.heuristic._score(game, move, player), move) for move in allowed
        ]
        best = max(score for score, _ in scored)
        return self.rng.choice(
            [move for score, move in scored if score >= best - 1e-9]
        )

    def _update_knowledge(self, game: Game, owner: Owner) -> None:
        if self.knowledge is None or self.knowledge.owner != owner:
            self.knowledge = OpponentKnowledge(owner)
            self.processed_records = 0
        # Set every call, not just on construction: a weight can differ between
        # two candidates sharing one agent class, and the flag has to follow it.
        self.knowledge.deduce_engineers = bool(
            self.heuristic.weights.use_engineer_deduction
        )
        for event in game.observations(owner, self.processed_records):
            self.knowledge.observe(event)
        self.processed_records = len(game.records)
        self.knowledge.forget_missing(
            {
                position
                for position, piece in game.board.items()
                if piece.owner == owner.other
            }
        )
        if self.heuristic.weights.use_rank_elimination:
            self.knowledge.eliminate_dead_ranks()

    def _acceptable(
        self, game: Game, move: Move, player: Owner, targets: list[Position]
    ) -> bool:
        target = game.board.get(move.dst)
        if target is None:
            return True
        if target.revealed or move.dst in targets:
            return True  # a flag or last-headquarters strike is always worth it
        attacker = game.board[move.src].kind
        assert self.knowledge is not None
        believed = self.knowledge.possible.get(move.dst)
        if believed:
            return _expected_battle(attacker, believed) > self.threshold
        if move.dst[0] in rear_rows(player.other):
            # Mine country, and only an engineer survives a mine.
            return attacker == PieceKind.ENGINEER
        return PRIOR_BATTLE[attacker] > self.threshold


def material_weights(base: BotWeights) -> BotWeights:
    """Pure trader: no headquarters plan, all value in captures."""
    return replace(
        base,
        hq_pressure=0.0,
        hq_strike=0.0,
        eval_hq_attack=0.0,
        eval_hq_attack_certain=0.0,
        eval_hq_defense=0.0,
        eval_hq_defense_certain=0.0,
        capture=base.capture * 2.0,
        eval_material=base.eval_material * 1.5,
    )


def defensive_weights(base: BotWeights) -> BotWeights:
    """Turtle: heavy own-flag defence, reluctant to trade forward."""
    return replace(
        base,
        forward=base.forward * 0.2,
        protect_flag=base.protect_flag * 6.0,
        eval_hq_defense=base.eval_hq_defense * 3.0,
        eval_hq_defense_certain=base.eval_hq_defense_certain * 3.0,
        hq_pressure=base.hq_pressure * 0.3,
    )


def engineer_preserving_weights(base: BotWeights) -> BotWeights:
    """Hoards engineers for mine clearing and refuses to spend them."""
    return replace(
        base,
        engineer_waste=base.engineer_waste * 3.0,
        engineer_mine=base.engineer_mine * 2.0,
        mine_risk=base.mine_risk * 2.0,
    )


@dataclass(frozen=True)
class AgentSpec:
    """A picklable recipe for one player."""

    name: str
    kind: str  # random | heuristic | search | hqrush | selective | oracle
    weights_path: str | None = None
    weights_style: str = "as_is"  # as_is | material | defensive | engineer
    samples: int = 3
    beam_width: int = 8
    reply_width: int = 4
    caution: float = 1.0
    noise: float | None = None
    #: Minimum expected battle outcome a `selective` agent will attack on.
    threshold: float = 0.05
    #: Fraction of hidden ranks an `oracle` agent sees. 1.0 is full
    #: clairvoyance; lower values are Suphx's dropout schedule.
    gamma: float = 1.0
    #: Fraction the oracle's *candidate generator* also sees. 0 keeps the
    #: beam honest, which is what makes the plain oracle a lower bound.
    gamma_beam: float = 0.0
    #: ``"module:attribute"`` naming a factory ``(spec, weights, seed) -> player``.
    #: Resolved by import *inside the worker*, so an agent that needs a network
    #: client or a third-party SDK can join the pool without `junqi` importing
    #: one -- this package is deliberately stdlib-only. A plain string keeps the
    #: spec picklable, which a callable would not be.
    builder: str | None = None
    #: Configuration for an external builder, as a tuple of pairs so the spec
    #: stays frozen and hashable.
    options: tuple[tuple[str, str], ...] = ()

    def build(self, weights: BotWeights, seed: int):
        if self.builder is not None:
            module_name, _, attribute = self.builder.partition(":")
            factory = getattr(import_module(module_name), attribute)
            return factory(self, weights, seed)
        if self.kind == "random":
            return RandomBot(seed=seed)
        if self.kind == "hqrush":
            return HQRushBot(seed=seed, caution=self.caution)
        shaped = STYLES[self.weights_style](weights)
        if self.noise is not None:
            shaped = replace(shaped, noise=self.noise)
        if self.kind == "selective":
            return SelectiveBot(shaped, seed=seed, threshold=self.threshold)
        if self.kind == "oracle":
            return OracleSearchBot(
                shaped,
                seed=seed,
                samples=self.samples,
                beam_width=self.beam_width,
                reply_width=self.reply_width,
                gamma=self.gamma,
                gamma_beam=self.gamma_beam,
            )
        if self.kind == "heuristic":
            return HeuristicBot(shaped, seed=seed)
        if self.kind == "search":
            return SearchBot(
                shaped,
                seed=seed,
                samples=self.samples,
                beam_width=self.beam_width,
                reply_width=self.reply_width,
            )
        raise ValueError(f"unknown agent kind: {self.kind}")


STYLES = {
    "as_is": lambda base: base,
    "material": material_weights,
    "defensive": defensive_weights,
    "engineer": engineer_preserving_weights,
}


@dataclass
class Pool:
    specs: list[AgentSpec] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.specs)


def standard_pool(
    stable: str | None = None,
    history: list[str] | None = None,
    anchor: str | None = None,
    exploiters: list[str] | None = None,
) -> Pool:
    """The default judging pool.

    `anchor` pins every weight-driven opponent to one fixed model file. Without
    it, `play_match` builds those opponents from *the subject's own weights*, so
    the opponents change with the model under test and the pool stops being a
    yardstick: a candidate would be measured against a distorted copy of itself
    while the baseline is measured against a distorted copy of *itself*. Always
    pass an anchor when comparing two models.
    """
    specs = [
        AgentSpec("random", "random"),
        AgentSpec("heuristic", "heuristic", weights_path=anchor),
        AgentSpec("heuristic-quiet", "heuristic", weights_path=anchor, noise=0.02),
        AgentSpec(
            "material", "heuristic", weights_path=anchor, weights_style="material"
        ),
        AgentSpec(
            "defensive", "heuristic", weights_path=anchor, weights_style="defensive"
        ),
        AgentSpec(
            "engineer", "heuristic", weights_path=anchor, weights_style="engineer"
        ),
        AgentSpec("hqrush", "hqrush"),
        AgentSpec("hqrush-careful", "hqrush", caution=2.0),
        # Attacks only on favourable odds. Without these the pool cannot tell a
        # well-priced attack term from a badly-priced one, because every other
        # weight-driven opponent carries the candidate's own attacking bias.
        AgentSpec("selective", "selective", weights_path=anchor),
        AgentSpec(
            "selective-strict", "selective", weights_path=anchor, threshold=0.35
        ),
        AgentSpec(
            "search-shallow",
            "search",
            weights_path=anchor,
            samples=1,
            beam_width=4,
            reply_width=1,
        ),
        AgentSpec(
            "search-mid",
            "search",
            weights_path=anchor,
            samples=3,
            beam_width=8,
            reply_width=3,
        ),
        AgentSpec(
            "search-deep",
            "search",
            weights_path=anchor,
            samples=6,
            beam_width=10,
            reply_width=4,
        ),
        # The pool's ceiling. Every other weight-driven member is the subject's
        # own policy class, so without this the pool cannot grade anything past
        # parity-with-self -- `search-mid` topped out near 60%. Cheating is the
        # point, and it is also cheap: at gamma 1 one sampled world suffices.
        AgentSpec("oracle", "oracle", weights_path=anchor, beam_width=10),
        AgentSpec(
            "oracle-half", "oracle", weights_path=anchor, samples=4,
            beam_width=10, gamma=0.5,
        ),
        # The genuine upper bound: the candidate generator cheats too, so a
        # move that is only good because the enemy ranks are visible can
        # actually reach the beam.
        AgentSpec(
            "oracle-perfect", "oracle", weights_path=anchor, beam_width=10,
            gamma=1.0, gamma_beam=1.0,
        ),
    ]
    if stable is not None:
        specs.append(
            AgentSpec("stable", "search", weights_path=stable, samples=3, beam_width=8)
        )
    # Best-response agents trained against one pool member (scripts/exploiter.py).
    # They are the other half of the ceiling problem: an oracle is stronger in a
    # way no human is, an exploiter is stronger in the way a *student of this bot*
    # would be.
    for index, path in enumerate(exploiters or []):
        specs.append(
            AgentSpec(
                f"exploiter-{Path(path).stem}",
                "search",
                weights_path=path,
                samples=3,
                beam_width=8,
            )
        )
    for index, path in enumerate(history or []):
        specs.append(
            AgentSpec(
                f"history-{index}",
                "search",
                weights_path=path,
                samples=3,
                beam_width=8,
            )
        )
    return Pool(specs)


def discover_exploiters(directory: str | Path, limit: int = 4) -> list[str]:
    """Accepted best-response models, newest first."""
    folder = Path(directory)
    if not folder.exists():
        return []
    found = sorted(
        folder.glob("*.json"), key=lambda path: path.stat().st_mtime, reverse=True
    )
    return [str(path) for path in found[:limit]]


def discover_history(directory: str | Path, limit: int = 3) -> list[str]:
    """Most recent archived models, newest first."""
    folder = Path(directory)
    if not folder.exists():
        return []
    snapshots = sorted(
        folder.glob("*.json"), key=lambda path: path.stat().st_mtime, reverse=True
    )
    return [str(path) for path in snapshots[:limit]]
