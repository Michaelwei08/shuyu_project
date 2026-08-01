from __future__ import annotations

import json
import random
from dataclasses import asdict, dataclass, fields
from pathlib import Path

from .board import CAMPS
from .deployment import rear_rows
from .game import Game, battle_outcome
from .types import Move, Owner, PieceKind, Position


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
    unknown_risk: float = 0.12
    belief_battle: float = 9.0
    hq_pressure: float = 1.2
    hq_strike: float = 14.0
    mine_risk: float = 0.5
    engineer_mine: float = 3.5
    engineer_waste: float = -4.5
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
    # A commander's death reveals its own flag, so its life buys concealment
    # that material value (11, against a general's 10) does not price in.
    # Tempting, and measurably wrong: at 18.0 this scored -0.0213 +/- 0.0136
    # over 814 paired games, i.e. worse. The arithmetic agrees -- a commander
    # attacking a non-rear square only loses to a bomb, about 8%, so the raid
    # is +9.6 in expectation and deterring it would need a coefficient near
    # 138, which would make the bot hoard its best piece. Left in at zero so
    # training can revisit it; do not hand-raise it without a paired result.
    eval_commander: float = 0.0

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

        targets = enemy_flag_squares(game, owner)
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
                score -= _piece_value(piece.kind) * weights.unknown_risk

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


def enemy_flag_squares(game: Game, owner: Owner) -> list[Position]:
    """Squares where ``owner``'s opponent may be hiding its flag.

    A revealed flag collapses this to one square; otherwise every still-held
    enemy headquarters is a candidate. Only occupancy is inspected, so this
    stays inside the hidden-rank invariant.
    """
    candidates = game.flag_candidates(owner.other)
    revealed = [square for square in candidates if game.board[square].revealed]
    return revealed or candidates


def _expected_battle(attacker: PieceKind, possible: frozenset[PieceKind]) -> float:
    if not possible:
        return 0.0
    return sum(battle_outcome(attacker, kind) for kind in possible) / len(possible)


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
