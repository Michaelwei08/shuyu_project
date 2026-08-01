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
from pathlib import Path

from .bot import BotWeights, HeuristicBot, _distance, _piece_value
from .game import Game
from .search_bot import SearchBot
from .types import Move, Owner, PieceKind


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
    kind: str  # "random" | "heuristic" | "search" | "hqrush"
    weights_path: str | None = None
    weights_style: str = "as_is"  # as_is | material | defensive | engineer
    samples: int = 3
    beam_width: int = 8
    reply_width: int = 4
    caution: float = 1.0
    noise: float | None = None

    def build(self, weights: BotWeights, seed: int):
        if self.kind == "random":
            return RandomBot(seed=seed)
        if self.kind == "hqrush":
            return HQRushBot(seed=seed, caution=self.caution)
        shaped = STYLES[self.weights_style](weights)
        if self.noise is not None:
            shaped = replace(shaped, noise=self.noise)
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
    ]
    if stable is not None:
        specs.append(
            AgentSpec("stable", "search", weights_path=stable, samples=3, beam_width=8)
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


def discover_history(directory: str | Path, limit: int = 3) -> list[str]:
    """Most recent archived models, newest first."""
    folder = Path(directory)
    if not folder.exists():
        return []
    snapshots = sorted(
        folder.glob("*.json"), key=lambda path: path.stat().st_mtime, reverse=True
    )
    return [str(path) for path in snapshots[:limit]]
