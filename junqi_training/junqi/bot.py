from __future__ import annotations

import json
import random
from dataclasses import asdict, dataclass, fields
from pathlib import Path

from .board import (
    CAMPS,
    HEADQUARTERS,
    engineer_rail_destinations,
    road_neighbors,
    straight_rail_destinations,
)
from .deployment import rear_rows
from .game import Game, battle_outcome
from .types import PIECE_COUNTS, Move, Owner, PieceKind, Position


@dataclass
class BotWeights:
    """The single source of truth for both engines' coefficients.

    Every tunable number in `junqi/` *and* in `web/lib/bot.ts` lives here.
    `junqi.web_export` writes `web/lib/weights.ts` from this dataclass, and
    `test_web_weights_are_in_sync` fails if the two drift apart.
    """

    # --- move heuristic -------------------------------------------------
    # There is deliberately no favourable/losing-battle pair here: `revealed`
    # is only ever set by the commander-death flag reveal, so a revealed piece
    # is always a flag and a "revealed non-flag battle" cannot occur. Those two
    # coefficients existed and were being mutated, costing the search two
    # dimensions of pure noise.
    capture: float = 2.8
    flag_capture: float = 120.0
    forward: float = 0.55
    camp: float = 1.1
    mobility: float = 0.08
    protect_flag: float = 0.3
    revealed_flag_hunt: float = 5.0
    belief_battle: float = 9.0
    # Superseded by `blind_battle`, and kept only so the paired harness can A/B
    # the two: `compare()` varies weights, not code, so deleting this branch
    # would apply the change to *both* sides of every comparison and cancel out.
    # Shipped at 0. Set it back to 0.0908 (and `blind_battle` to 0) to recover
    # the pre-2026-08-01 behaviour exactly.
    unknown_risk: float = 0.0
    hq_pressure: float = 1.2
    hq_strike: float = 14.0
    mine_risk: float = 0.5
    engineer_mine: float = 3.5
    # An engineer beats a mine, ties a bomb, and loses to all nine other
    # ranks. Attacking anything that is not plausibly a mine is simply handing
    # the piece over -- and it is one of only three answers to a mine. At -4.5
    # this merely cancelled the capture bonus (+3.2), leaving the bot
    # indifferent to throwing engineers away.
    engineer_waste: float = -12.0
    # Attacking a square we know nothing about. This used to be
    # `-piece_value * unknown_risk`, which charged for *what the attacker is
    # worth* -- but what decides a blind attack is *how likely it is to lose*,
    # and those run in opposite directions. The old term therefore scored the
    # commander lowest (+1.80) and the engineer highest (+2.53) for the same
    # blind attack, while ten replayed games had the commander at 3W/1T/0L and
    # the engineer at 1W/2T/5L. It paid the bot to probe with the pieces that
    # lose. Now the branch prices the expected battle outcome against the prior
    # over plausible surviving ranks, exactly as `belief_battle` does for a
    # square we *have* deduced something about -- separate coefficient because a
    # deduced set is far better information than a prior and the two should not
    # be forced equal.
    blind_battle: float = 9.0
    noise: float = 0.18

    # --- state evaluation -----------------------------------------------
    eval_material: float = 1.7
    eval_mobility: float = 0.06
    eval_terminal: float = 2_000.0
    eval_hq_attack: float = 0.9
    eval_hq_attack_certain: float = 2.6
    eval_hq_defense: float = 1.1
    eval_hq_defense_certain: float = 3.0
    # The distance terms above are linear over 0..12, so "they take my flag next
    # ply" scored only 1.8x "they are six squares away" -- less than a single
    # capture. These two are the sharp part of the signal.
    eval_hq_breach: float = 26.0
    eval_hq_guard: float = 5.5
    # Replayed losses were not the bot ignoring a breach -- in nine of nine it
    # answered every breach it *could*, and on the fatal ply it had no legal
    # move onto the approach square at all. `eval_hq_breach` cannot help there:
    # with no answer available the penalty is identical across every candidate
    # and cancels out. Worse, a search that assumes it answers next ply prices
    # spending its last defender as free. This term is about *supply* rather
    # than the current occupant, so a continuation cannot dissolve it: while an
    # enemy is within two moves of a live headquarters, having fewer than two
    # pieces covering its approaches is charged for.
    eval_hq_supply: float = 0.0
    # How much the fast move heuristic counts next to the sampled rollout when
    # `SearchBot` picks a move. The two engines silently disagreed on this:
    # Python used `base_score + rollout` while `web/lib/bot.ts` used
    # `base * 2 + search`, so identical weights produced a different policy in
    # the browser than in the harness that measured them. Shipped at 1.0, which
    # keeps the Python side -- the side every result was measured on -- exactly
    # as it was, and brings the browser into line. Which value is actually
    # better is now a weight the harness can answer.
    search_base_weight: float = 1.0
    # How much the *modelled* opponent inside a rollout knows about our army.
    # It knew nothing: `SearchBot._rollout` builds a `HeuristicBot` for reply
    # ranking and never assigns `.knowledge`, so every sampled world assumed an
    # opponent with no deductions at all while the bot itself runs a full
    # `OpponentKnowledge`. The search therefore cannot represent an opponent who
    # has worked out what our pieces are, and systematically over-rates lines
    # that only survive against a blind reply.
    #
    # 0 reproduces that exactly; 1 models an omniscient opponent. Neither
    # extreme is obviously right -- the oracle diagnostic puts perfect
    # information at about +10 points of win rate, so assuming an omniscient
    # replier overstates the danger by roughly that much. Ships at 0 so the
    # harness decides. Note this reads only *our own* ranks, never theirs.
    reply_insight: float = 0.0
    # Weight on `models/value.json`, the leaf evaluation fit to self-play
    # outcomes rather than tuned by hand (see `junqi/value.py`). The learned
    # model is *added* to the existing terms rather than replacing them, so this
    # sweeps continuously from "exactly today's bot" at 0 to "the learned model
    # dominates" -- which keeps the whole change inside what `compare()` can
    # measure as a weight. Ships at 0; nothing is adopted before a paired run
    # says so. Inert with no `models/value.json` on disk.
    eval_value_scale: float = 0.0
    # A commander's death reveals its own flag, so its life buys concealment
    # that material value (11, against a general's 10) does not price in.
    # Tempting, and measurably wrong: at 18.0 this scored -0.0213 +/- 0.0136
    # over 814 paired games, i.e. worse. The arithmetic agrees -- a commander
    # attacking a non-rear square only loses to a bomb, about 8%, so the raid
    # is +9.6 in expectation and deterring it would need a coefficient near
    # 138, which would make the bot hoard its best piece. Left in at zero so
    # training can revisit it; do not hand-raise it without a paired result.
    eval_commander: float = 0.0
    # Only engineers turn corners on the railway, so making one of those turns
    # announces the piece's rank. Spending that disguise for nothing walks the
    # engineer into an easy capture -- and engineers are the only answer to a
    # mine, so they are not spare material.
    engineer_expose: float = -6.0
    # Taking the flag is not the only way to win: a side with no legal move
    # loses. Mines and flags never move and a headquarters piece is frozen, so
    # an opponent can be reduced to zero *mobile* pieces while still holding
    # six. Material alone barely notices; this makes the last few captures
    # worth what they actually are.
    eval_immobilize: float = 110.0
    # 1 = count threat in moves (railways make E2 two moves from B1 while
    # Manhattan calls it four); 0 = the old Manhattan metric. A weight rather
    # than a constant purely so the paired harness can A/B it: a code-level
    # switch would apply to both sides of a comparison and cancel out.
    use_move_distance: float = 1.0
    # Where the headquarters distance terms fall to zero. Manhattan spans 0..15
    # so 12 suited it; move distance spans 0..5, which is why this was dropped
    # to 6. Both went in together and may pull opposite ways, so it is a weight
    # too -- the metric and its rescaling need separating.
    eval_horizon: float = 6.0
    # Cost of ending a move on a headquarters square, per point of piece value.
    # Immobility attaches to the *square*, not to the piece deployed there --
    # `_destinations` returns an empty set for any source in HEADQUARTERS -- so
    # a piece that wins a headquarters probe is alive, still counted in
    # material, and permanently useless.
    #
    # Nothing priced that before. The mobility term is the only one that could
    # notice, and it is skipped whenever the target is an unrevealed occupied
    # square, which is exactly what a probe is; `_rollout` runs `quick=True`,
    # which drops it entirely. Meanwhile `blind_battle` correctly prefers the
    # attacker with the best odds. The two compose into "send the commander",
    # because *winning* the probe is what entombs you -- so the stronger the
    # prober, the likelier the worst non-winning outcome.
    #
    # Measured over 320 pool games: 136 entombments, mean value 8.8, and
    # COMMANDER was the single most frequent rank to be entombed (41 of 136)
    # against ENGINEER's 1. Scaling by `_piece_value` is right here for the
    # reason it was wrong in `blind_battle`: this term prices lost future value,
    # not odds.
    #
    # Ships at 0 until it has a paired p-value. `models/ab/entomb-on.json` is
    # the switched-on end.
    eval_hq_entomb: float = 0.0
    # 1 = also run `OpponentKnowledge.eliminate_dead_ranks`, which counts the
    # casualty record and drops ranks it proves extinct from every belief set;
    # 0 = per-square battle constraints only, the deduction as it has always
    # been. A weight for the same reason as `use_move_distance`: `compare()`
    # varies weights, not code, so a plain `if` would apply to both sides and
    # cancel. Ships at 0 until it has a paired p-value (D027).
    use_rank_elimination: float = 0.0

    @classmethod
    def load(cls, path: str | Path) -> "BotWeights":
        with Path(path).open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
        known = {descriptor.name for descriptor in fields(cls)}
        # Tolerate models written before a weight was added or removed.
        return cls(**{key: value for key, value in payload.items() if key in known})

    def save(self, path: str | Path) -> None:
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("w", encoding="utf-8") as handle:
            json.dump(asdict(self), handle, ensure_ascii=False, indent=2)
            handle.write("\n")


class HeuristicBot:
    def __init__(
        self, weights: BotWeights | None = None, seed: int | None = None
    ) -> None:
        self.weights = weights or BotWeights()
        self.rng = random.Random(seed)
        # Set by SearchBot so move ordering can use privately deduced ranks.
        self.knowledge: dict[Position, frozenset[PieceKind]] | None = None

    def choose_move(self, game: Game, owner: Owner | None = None) -> Move:
        player = game.turn if owner is None else owner
        moves = game.legal_moves(player)
        if not moves:
            raise ValueError("当前玩家没有合法走法")
        scored = [(self._score(game, move, player), move) for move in moves]
        best = max(score for score, _ in scored)
        candidates = [move for score, move in scored if score >= best - 1e-9]
        return self.rng.choice(candidates)

    def _score(
        self, game: Game, move: Move, owner: Owner, quick: bool = False
    ) -> float:
        """Score one move.

        `quick` drops the mobility term, which is the only part that has to
        simulate the move. That simulation dominated search time when used to
        rank the opponent's replies inside a rollout, where all we need is a
        rough top-K ordering.
        """
        piece = game.board[move.src]
        target = game.board.get(move.dst)
        weights = self.weights
        score = self.rng.uniform(-weights.noise, weights.noise)

        direction = 1 if owner == Owner.BOT else -1
        score += direction * (move.dst[0] - move.src[0]) * weights.forward
        if move.dst in CAMPS:
            score += weights.camp
        if piece.kind == PieceKind.ENGINEER and _reveals_engineer(game, move):
            score += weights.engineer_expose

        targets = enemy_flag_squares(game, owner, self.knowledge)
        certain = len(targets) == 1

        if target is not None:
            score += weights.capture
            believed = None if self.knowledge is None else self.knowledge.get(move.dst)
            if target.revealed or move.dst in targets:
                # `revealed` means a flag; a still-held enemy headquarters holds
                # either the flag or the last decoy hiding it.
                score += (
                    weights.flag_capture
                    if target.revealed or certain
                    else weights.hq_strike
                )
            elif believed:
                score += _expected_battle(piece.kind, believed) * weights.belief_battle
            elif move.dst[0] in rear_rows(owner.other):
                # Mines only ever sit in the rear two rows, and only an
                # engineer survives one.
                score += (
                    weights.engineer_mine
                    if piece.kind == PieceKind.ENGINEER
                    else -_piece_value(piece.kind) * weights.mine_risk
                )
            elif piece.kind == PieceKind.ENGINEER:
                # An engineer loses to every rank except a mine.
                score += weights.engineer_waste
            else:
                score += _prior_battle(piece.kind) * weights.blind_battle
                score -= _piece_value(piece.kind) * weights.unknown_risk

        # Ending on a headquarters freezes the piece for the rest of the game.
        # Skipped when the square is known to hold the flag, because then the
        # game ends on this move and being frozen costs nothing. Cheap enough
        # for `quick=True`, unlike the mobility term that should have caught it.
        takes_flag = target is not None and (
            target.revealed or (certain and move.dst in targets)
        )
        if weights.eval_hq_entomb and move.dst in HEADQUARTERS and not takes_flag:
            score -= _piece_value(piece.kind) * weights.eval_hq_entomb

        if not quick and weights.mobility and (target is None or target.revealed):
            simulated = game.clone()
            simulated.turn = owner
            simulated.apply(move)
            if not simulated.over:
                score += len(simulated.legal_moves(owner)) * weights.mobility

        own_flag = game.flag_position(owner)
        if own_flag is not None:
            before = _distance(move.src, own_flag)
            after = _distance(move.dst, own_flag)
            score += (before - after) * weights.protect_flag
        if targets:
            before = min(_distance(move.src, square) for square in targets)
            after = min(_distance(move.dst, square) for square in targets)
            hunt = weights.revealed_flag_hunt if certain else weights.hq_pressure
            score += (before - after) * hunt
        return score


def enemy_flag_squares(
    game: Game,
    owner: Owner,
    knowledge: dict[Position, frozenset[PieceKind]] | None = None,
) -> list[Position]:
    """Squares where ``owner``'s opponent may be hiding its flag.

    A revealed flag collapses this to one square. Failing that, a *failed*
    attack on a headquarters proves the square is not the flag -- a flag loses
    to every rank -- so the other headquarters holds it with certainty. That
    deduction is already in ``knowledge``: a lost attack records the ranks that
    beat our piece, and FLAG beats nothing, so it is absent from that set.
    Occupancy and our own battle history only, so no hidden rank is read.
    """
    candidates = game.flag_candidates(owner.other)
    revealed = [square for square in candidates if game.board[square].revealed]
    if revealed:
        return revealed
    if knowledge:
        possible = [
            square
            for square in candidates
            if PieceKind.FLAG in knowledge.get(square, frozenset({PieceKind.FLAG}))
        ]
        if possible:
            return possible
    return candidates


def _reveals_engineer(game: Game, move: Move) -> bool:
    """True when only an engineer could have made this move.

    Engineers alone turn corners on the railway, so such a move tells the
    opponent exactly what the piece is.
    """
    if move.dst in road_neighbors(move.src):
        return False
    occupied = set(game.board)
    if move.dst in straight_rail_destinations(move.src, occupied):
        return False
    return move.dst in engineer_rail_destinations(move.src, occupied)


def _expected_battle(attacker: PieceKind, possible: frozenset[PieceKind]) -> float:
    if not possible:
        return 0.0
    return sum(battle_outcome(attacker, kind) for kind in possible) / len(possible)


def _build_prior_battle() -> dict[PieceKind, float]:
    """Expected battle outcome against a square we know nothing about.

    A flag never leaves a headquarters and a mine never leaves the rear two
    rows, so neither can be the occupant of the square this branch scores --
    the prior is the rest of the army, weighted by how many of each rank exist.
    Static, so build it once at import rather than per scored move.
    """
    pool = {
        kind: count
        for kind, count in PIECE_COUNTS.items()
        if kind not in (PieceKind.FLAG, PieceKind.MINE)
    }
    total = sum(pool.values())
    return {
        attacker: sum(
            battle_outcome(attacker, defender) * count
            for defender, count in pool.items()
        )
        / total
        for attacker in PieceKind
    }


#: Ranges from +0.86 (commander) to -0.76 (engineer): monotone in rank, which is
#: the property the old `unknown_risk` term got backwards.
PRIOR_BATTLE: dict[PieceKind, float] = _build_prior_battle()


def _prior_battle(attacker: PieceKind) -> float:
    return PRIOR_BATTLE[attacker]


def _distance(left: tuple[int, int], right: tuple[int, int]) -> int:
    return abs(left[0] - right[0]) + abs(left[1] - right[1])


def _piece_value(kind: PieceKind) -> float:
    if kind == PieceKind.FLAG:
        return 50.0
    if kind == PieceKind.BOMB:
        return 7.0
    if kind == PieceKind.MINE:
        return 5.0
    return float(12 - kind.value)
