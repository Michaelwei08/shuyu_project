"""Fit `models/value.json` to self-play outcomes.

    python -m junqi.value_training --games 600 --workers 20

Three things here are load-bearing, and each of them is a way the naive version
of this gets a flattering number that does not survive contact with the search.

**Features are computed on a determinized world, not the true board.**
`SearchBot._state_value` only ever sees sampled worlds, where the opponent's
ranks are guesses. Training on the true board would fit "how many engineers do
they actually have left", which is not a question the deployed evaluation can
ask. Own composition, all positions, mobility and headquarters occupancy are
exact either way; only the opponent's rank composition is sampled.

**The train/test split is by game, never by position.** Plies within one game
share an outcome and are hugely autocorrelated, so a random row split leaks the
label and inflates AUC by roughly 0.15 -- the same grouped-split lesson the
Raman project learned the hard way with `fold_masks`.

**Draws are dropped rather than labelled 0.5.** A draw carries no information
about which side was better and labelling it at the midpoint just shrinks every
coefficient toward zero.
"""

from __future__ import annotations

import argparse
import random
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from .arena import MAX_MOVES, default_workers, make_opening
from .bot import BotWeights
from .game import Game
from .opponents import standard_pool
from .search_bot import SearchBot
from .types import Owner
from .value import auc, features, fit, phase_of

DEFAULT_MODEL = Path("models/bot_weights.json")
DEFAULT_OUTPUT = Path("models/value.json")
DEFAULT_ANCHOR = Path("models/defaults.json")


Row = tuple[list[float], float, int, float]


def _collect(task: tuple[int, str, str | None, frozenset[str], int]) -> list[Row]:
    """Play one game and return `(features, label, phase)` for every ply.

    Labels are from the *bot* side's point of view; both sides are recorded, the
    human side with its own perspective and its own outcome, so one game yields
    two views and the model never learns "north wins".
    """
    seed, model_path, anchor, excluded, stride = task
    weights = BotWeights.load(model_path)
    specs = [
        spec for spec in standard_pool(anchor=anchor).specs
        if spec.name not in excluded
    ]
    if not specs:
        raise SystemExit("every pool opponent was excluded")
    spec = specs[seed % len(specs)]

    game = Game(board=make_opening(seed), turn=Owner(seed % 2))
    subject = SearchBot(weights, seed=seed * 2 + 1, samples=3, beam_width=8,
                        reply_width=4)
    opponent_weights = (
        BotWeights.load(spec.weights_path) if spec.weights_path else weights
    )
    challenger = spec.build(opponent_weights, seed=seed * 2)
    players = {Owner.BOT: subject, Owner.HUMAN: challenger}

    rows: list[tuple[list[float], float, int, Owner]] = []  # side kept for labelling
    sampler = SearchBot(weights, seed=seed * 3 + 7, samples=1, beam_width=1)
    rng = random.Random(seed * 7919 + 13)
    while not game.over and game.move_count < MAX_MOVES:
        # Subsample plies. Consecutive positions in one game are nearly
        # identical and share a label, so most of them are redundant; this
        # is what keeps a 6400-game collection inside memory, and it weakens
        # the within-game correlation into the bargain.
        for side in Owner if game.move_count % stride == 0 else ():
            # Exactly what the search sees: the opponent's hidden ranks resampled.
            world = sampler._determinize(game, side.other, rng)
            rows.append(
                (
                    features(world, side),
                    # The incumbent scored on the *same* world, so the
                    # before/after is a like-for-like comparison rather than the
                    # learned model being handed better inputs.
                    sampler._state_value(world, side),
                    phase_of(game.move_count),
                    side,
                )
            )
        game.apply(players[game.turn].choose_move(game))

    if game.winner is None:
        return []  # a draw says nothing about who stood better
    return [
        (vector, 1.0 if side == game.winner else 0.0, phase, shipped)
        for vector, shipped, phase, side in rows
    ]


def collect(
    games: int,
    model: Path,
    anchor: Path | None,
    workers: int | None,
    excluded: frozenset[str] = frozenset(),
    stride: int = 1,
    offset: int = 900_001,
) -> list[list[Row]]:
    tasks = [
        (seed, str(model), str(anchor) if anchor else None, excluded, stride)
        for seed in range(offset, offset + games)
    ]
    count = default_workers() if workers is None else workers
    if count <= 1:
        return [_collect(task) for task in tasks]
    with ProcessPoolExecutor(max_workers=count) as executor:
        return list(executor.map(_collect, tasks, chunksize=2))


def report(name: str, scores: list[float], rows: list[Row]) -> str:
    line = f"{name:<34}"
    for phase in range(3):
        subset = [
            (score, label)
            for score, (_, label, bucket, _shipped) in zip(scores, rows, strict=True)
            if bucket == phase
        ]
        if len(subset) < 200:
            line += f"{'--':>10}"
            continue
        line += f"{auc([s for s, _ in subset], [y for _, y in subset]):>10.3f}"
    line += f"{auc(scores, [y for _, y, _, _ in rows]):>10.3f}"
    return line


def main() -> None:
    parser = argparse.ArgumentParser(description="用自对弈结果拟合局面价值函数")
    parser.add_argument("--games", type=int, default=600)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--anchor", type=Path, default=None)
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--epochs", type=int, default=40)
    parser.add_argument(
        "--stride",
        type=int,
        default=1,
        help="每隔几手取一个训练样本。相邻局面几乎一样且共享标签，大规模收集时用 4 或 8 可以大幅降内存，信息损失很小",
    )
    parser.add_argument(
        "--holdout", type=float, default=0.3, help="按对局（不是按局面）留出的比例"
    )
    parser.add_argument(
        "--exclude",
        nargs="*",
        default=[],
        metavar="OPPONENT",
        help=(
            "训练时不与这些对手对局。这是**对手层面**的留出，不是局面层面的。"
            "价值函数拟合的是“对着这个对手池，棋局最后怎么收场”，"
            "再拿同一个池子评测，涨幅里有多少是棋力、有多少是过拟合这个池子，"
            "分不开。先 --exclude 掉几个对手训练、再用它们评测，才是干净的检验"
        ),
    )
    arguments = parser.parse_args()
    anchor = arguments.anchor
    if anchor is None and DEFAULT_ANCHOR.exists():
        anchor = DEFAULT_ANCHOR

    print(f"collecting {arguments.games} games ...", flush=True)
    excluded = frozenset(arguments.exclude)
    if excluded:
        print(f"holding out opponents: {sorted(excluded)}")
    games = collect(
        arguments.games,
        arguments.model,
        anchor,
        arguments.workers,
        excluded,
        arguments.stride,
    )
    games = [rows for rows in games if rows]
    positions = sum(len(rows) for rows in games)
    print(f"{len(games)} decisive games, {positions} labelled positions")
    if not games:
        raise SystemExit("no decisive games collected")

    # Split by GAME. A row split leaks the outcome across plies of one game.
    order = list(range(len(games)))
    random.Random(11).shuffle(order)
    cut = int(len(order) * (1.0 - arguments.holdout))
    train = [row for index in order[:cut] for row in games[index]]
    test = [row for index in order[cut:] for row in games[index]]
    print(f"{cut} train / {len(order) - cut} test games; "
          f"{len(train)} / {len(test)} positions")

    model = fit([(v, y) for v, y, _, _ in train], epochs=arguments.epochs)

    print(f"\n{'AUC by phase':<34}{'0-20':>10}{'20-45':>10}{'45+':>10}{'all':>10}")
    print("-" * 74)
    # Both scored on the *same* held-out worlds, so this is like for like: the
    # learned model is not being handed inputs the incumbent never saw.
    print(report("shipped _state_value", [s for _, _, _, s in test], test))
    print(
        report(
            "learned value model",
            [model.predict_features(v) for v, _, _, _ in test],
            test,
        )
    )

    model.save(arguments.output)
    print(f"\nsaved {arguments.output}")
    print("next: python -m junqi.web_export, then benchmark with "
          "--baseline models/ab/pre-value.json")


if __name__ == "__main__":
    main()
