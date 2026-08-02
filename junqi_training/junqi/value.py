"""A leaf evaluation fit to self-play outcomes instead of tuned by hand.

Four rounds of coefficient search, roughly 30,000 paired games, all came back
within noise of zero. An oracle diagnostic then priced the three subsystems:
perfect knowledge of every hidden rank is worth about +10 points of win rate,
quadrupling the search budget about +5, and coefficients about nothing. What
the coefficients *describe* -- the leaf evaluation -- had never been fit to
anything, and it spends nine of its terms on headquarters distance while having
no term at all for army composition.

So: same shape as `BotWeights` (a flat list of floats, one JSON file, exported
to the browser), but the numbers come from logistic regression against the
eventual result rather than from an evolutionary search over win/loss bits. The
~30,000 games already played contained millions of labelled positions and threw
all of them away.

The model is deliberately linear and small. It has to run in a browser inside a
420ms move budget with no dependencies, and a dot product is cheaper than the
`legalMoves`-twice evaluation it sits beside.

Every feature is replicated once per phase bucket, with only the active bucket
non-zero -- that is an interaction with game phase at the cost of a longer,
sparser dot product. It is not decoration: composition features measured
*worse* than useless in the opening and dominant after ply 45, so one coefficient
per feature cannot fit both.
"""

from __future__ import annotations

import json
import math
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

from .board import HEADQUARTERS, move_distance
from .game import Game
from .types import Owner, PieceKind, Position

#: Ply boundaries between phase buckets. Games run ~80 plies; the split is where
#: the measured feature importances change character rather than anywhere
#: principled.
PHASE_EDGES = (20, 45)
PHASES = len(PHASE_EDGES) + 1

#: Ranks worth counting individually. The rest carry little signal beyond their
#: contribution to material, and every extra feature is a parameter to fit.
TRACKED: tuple[PieceKind, ...] = (
    PieceKind.COMMANDER,
    PieceKind.GENERAL,
    PieceKind.ENGINEER,
    PieceKind.BOMB,
    PieceKind.MINE,
)

#: Where the raider-distance one-hot stops caring.
REACH_BUCKETS = 6

_VALUES: dict[PieceKind, float] = {
    PieceKind.FLAG: 50.0,
    PieceKind.BOMB: 7.0,
    PieceKind.MINE: 5.0,
}


def _value(kind: PieceKind) -> float:
    return _VALUES.get(kind, float(12 - kind.value))


def phase_of(move_count: int) -> int:
    for index, edge in enumerate(PHASE_EDGES):
        if move_count < edge:
            return index
    return PHASES - 1


def _mobile(game: Game, side: Owner) -> int:
    return sum(
        1
        for position, piece in game.board.items()
        if piece.owner == side
        and piece.kind.movable
        and position not in HEADQUARTERS
    )


def _across_river(game: Game, side: Owner) -> int:
    """Pieces that have committed to the far half of the board."""
    return sum(
        1
        for position, piece in game.board.items()
        if piece.owner == side
        and (position[0] >= 6 if side == Owner.BOT else position[0] <= 5)
    )


def _reach(game: Game, side: Owner, targets: list[Position]) -> int:
    if not targets:
        return REACH_BUCKETS - 1
    best = min(
        (
            move_distance(position, target)
            for position, piece in game.board.items()
            if piece.owner == side
            and piece.kind.movable
            and position not in HEADQUARTERS
            for target in targets
        ),
        default=REACH_BUCKETS - 1,
    )
    return min(best, REACH_BUCKETS - 1)


def base_features(game: Game, owner: Owner) -> list[float]:
    """Features from ``owner``'s point of view.

    Read off whatever board is passed in. `SearchBot` evaluates *determinized*
    worlds, so at inference the opponent's ranks are sampled rather than known
    -- which means training has to sample them the same way or the composition
    features mean something different in the two settings. `value_training`
    does exactly that.
    """
    enemy = owner.other
    own_pieces = [piece for piece in game.board.values() if piece.owner == owner]
    enemy_pieces = [piece for piece in game.board.values() if piece.owner == enemy]
    own_counts: dict[PieceKind, int] = {}
    enemy_counts: dict[PieceKind, int] = {}
    for piece in own_pieces:
        own_counts[piece.kind] = own_counts.get(piece.kind, 0) + 1
    for piece in enemy_pieces:
        enemy_counts[piece.kind] = enemy_counts.get(piece.kind, 0) + 1

    material = sum(_value(piece.kind) for piece in own_pieces) - sum(
        _value(piece.kind) for piece in enemy_pieces
    )
    own_hq = game.flag_candidates(owner)
    enemy_hq = game.flag_candidates(enemy)

    features = [
        1.0,
        material / 120.0,
        len(own_pieces) / 25.0,
        len(enemy_pieces) / 25.0,
        _mobile(game, owner) / 25.0,
        _mobile(game, enemy) / 25.0,
        len(own_hq) / 2.0,
        len(enemy_hq) / 2.0,
        _across_river(game, owner) / 25.0,
        _across_river(game, enemy) / 25.0,
    ]
    for kind in TRACKED:
        features.append(own_counts.get(kind, 0) / 3.0)
        features.append(enemy_counts.get(kind, 0) / 3.0)

    # One-hot rather than linear: the measured effect of a raider closing on a
    # headquarters is a step at "adjacent", not a slope.
    attack = _reach(game, owner, enemy_hq)
    defend = _reach(game, enemy, own_hq)
    features.extend(1.0 if index == attack else 0.0 for index in range(REACH_BUCKETS))
    features.extend(1.0 if index == defend else 0.0 for index in range(REACH_BUCKETS))
    return features


BASE_WIDTH = len(base_features(Game(board={}), Owner.BOT))
WIDTH = BASE_WIDTH * PHASES


def features(game: Game, owner: Owner) -> list[float]:
    """The full sparse vector: base features, in the active phase's slot."""
    vector = [0.0] * WIDTH
    offset = phase_of(game.move_count) * BASE_WIDTH
    for index, value in enumerate(base_features(game, owner)):
        vector[offset + index] = value
    return vector


@dataclass
class ValueModel:
    """P(owner eventually wins), as a logistic function of `features`."""

    weights: list[float] = field(default_factory=lambda: [0.0] * WIDTH)

    def predict_features(self, vector: list[float]) -> float:
        total = sum(w * x for w, x in zip(self.weights, vector, strict=True) if x)
        if total >= 0.0:
            return 1.0 / (1.0 + math.exp(-min(total, 60.0)))
        return math.exp(max(total, -60.0)) / (1.0 + math.exp(max(total, -60.0)))

    def predict(self, game: Game, owner: Owner) -> float:
        return self.predict_features(features(game, owner))

    def advantage(self, game: Game, owner: Owner) -> float:
        """Centred on 0 so the term is symmetric and adds cleanly to the rest."""
        return 2.0 * self.predict(game, owner) - 1.0

    @classmethod
    def load(cls, path: str | Path) -> "ValueModel":
        with Path(path).open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
        weights = list(payload["weights"])
        if len(weights) != WIDTH:
            raise ValueError(
                f"value model has {len(weights)} weights, this build expects "
                f"{WIDTH} -- refit with `python -m junqi.value_training`"
            )
        return cls(weights=weights)

    def save(self, path: str | Path) -> None:
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "base_width": BASE_WIDTH,
            "phases": PHASES,
            "phase_edges": list(PHASE_EDGES),
            "weights": self.weights,
        }
        with destination.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)
            handle.write("\n")


DEFAULT_PATH = Path("models/value.json")
_CACHE: dict[str, ValueModel | None] = {}


def load_default(path: str | Path = DEFAULT_PATH) -> ValueModel | None:
    """The shipped model, or ``None`` when there is not one yet.

    Cached per process because every worker in an evaluation run loads the same
    file, and `_state_value` is on the hottest path in the search.
    """
    key = str(path)
    if key not in _CACHE:
        candidate = Path(path)
        try:
            _CACHE[key] = ValueModel.load(candidate) if candidate.exists() else None
        except (ValueError, KeyError, json.JSONDecodeError):
            _CACHE[key] = None
    return _CACHE[key]


def fit(
    rows: list[tuple[list[float], float]],
    epochs: int = 40,
    rate: float = 0.5,
    decay: float = 1e-6,
    seed: int = 0,
) -> ValueModel:
    """Plain SGD logistic regression. Stdlib only, no third-party dependency.

    `rows` are `(feature vector, label)` with the label 1 for an eventual win
    and 0 for a loss. Draws are dropped by the caller rather than labelled 0.5,
    which would pull every coefficient toward zero for no information.
    """
    import random as _random

    rng = _random.Random(seed)
    model = ValueModel()
    order = list(range(len(rows)))
    for epoch in range(epochs):
        rng.shuffle(order)
        step = rate / (1.0 + epoch)
        for index in order:
            vector, label = rows[index]
            error = model.predict_features(vector) - label
            for position, value in enumerate(vector):
                if value:
                    model.weights[position] -= step * (
                        error * value + decay * model.weights[position]
                    )
    return model


def auc(scores: list[float], labels: list[float]) -> float:
    """Probability a random win outranks a random loss. Ties count a half."""
    positives = [s for s, y in zip(scores, labels, strict=True) if y > 0.5]
    negatives = [s for s, y in zip(scores, labels, strict=True) if y <= 0.5]
    if not positives or not negatives:
        return 0.5
    ordered = sorted(
        [(s, 1) for s in positives] + [(s, 0) for s in negatives], key=lambda r: r[0]
    )
    rank_sum = 0.0
    index = 0
    while index < len(ordered):
        stop = index
        while stop + 1 < len(ordered) and ordered[stop + 1][0] == ordered[index][0]:
            stop += 1
        average = (index + stop) / 2.0 + 1.0
        for position in range(index, stop + 1):
            if ordered[position][1]:
                rank_sum += average
        index = stop + 1
    count_p, count_n = len(positives), len(negatives)
    return (rank_sum - count_p * (count_p + 1) / 2.0) / (count_p * count_n)


FeatureFn = Callable[[Game, Owner], list[float]]
