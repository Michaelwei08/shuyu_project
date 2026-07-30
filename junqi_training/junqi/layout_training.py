from __future__ import annotations

import argparse
import random
from pathlib import Path

from .bot import BotWeights, HeuristicBot
from .deployment import (
    random_deployment,
    save_deployment,
    swap_pieces,
)
from .game import Game
from .training import DEFAULT_MODEL
from .types import Owner, Piece, Position

DEFAULT_LAYOUT = Path("models/bot_deployment.json")


def play_layout(
    layout: dict[Position, Piece],
    weights: BotWeights,
    seed: int,
    max_moves: int = 240,
) -> float:
    game = Game.new(
        seed=seed,
        first=Owner(seed % 2),
        bot_deployment=layout,
    )
    players = {
        Owner.HUMAN: HeuristicBot(weights, seed=seed * 2 + 1),
        Owner.BOT: HeuristicBot(weights, seed=seed * 2 + 2),
    }
    while not game.over and game.move_count < max_moves:
        game.apply(players[game.turn].choose_move(game))
    if game.winner == Owner.BOT:
        return 1.0
    if game.winner == Owner.HUMAN:
        return 0.0
    return 0.5


def score_layout(
    layout: dict[Position, Piece],
    weights: BotWeights,
    seeds: list[int],
) -> float:
    return sum(play_layout(layout, weights, seed) for seed in seeds) / len(seeds)


def mutate_layout(
    layout: dict[Position, Piece], rng: random.Random, swaps: int = 2
) -> dict[Position, Piece]:
    candidate = layout.copy()
    positions = list(candidate)
    completed = 0
    attempts = 0
    while completed < swaps and attempts < 100:
        attempts += 1
        left, right = rng.sample(positions, 2)
        try:
            swap_pieces(candidate, Owner.BOT, left, right)
        except ValueError:
            continue
        completed += 1
    if completed < swaps:
        raise RuntimeError("无法生成合法布局变体")
    return candidate


def train_layout(
    generations: int,
    games: int,
    seed: int,
    weights: BotWeights,
    output: Path,
) -> dict[Position, Piece]:
    rng = random.Random(seed)
    incumbent = random_deployment(Owner.BOT, rng)
    for generation in range(1, generations + 1):
        candidate = mutate_layout(incumbent, rng)
        seeds = [seed + generation * 10_000 + index for index in range(games)]
        incumbent_score = score_layout(incumbent, weights, seeds)
        candidate_score = score_layout(candidate, weights, seeds)
        accepted = candidate_score > incumbent_score
        if accepted:
            incumbent = candidate
        print(
            f"布局 {generation:02d}/{generations} | "
            f"原布局 {incumbent_score:.1%} | 候选 {candidate_score:.1%} | "
            f"{'采用' if accepted else '保留'}"
        )
    save_deployment(incumbent, output)
    print(f"布局模型已保存：{output.resolve()}")
    return incumbent


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="训练 Bot 初始布局")
    parser.add_argument("--generations", type=int, default=10)
    parser.add_argument("--games", type=int, default=6)
    parser.add_argument("--seed", type=int, default=20260802)
    parser.add_argument("--weights", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--output", type=Path, default=DEFAULT_LAYOUT)
    return parser


def main() -> None:
    arguments = build_parser().parse_args()
    if arguments.generations < 1 or arguments.games < 2:
        raise SystemExit("generations 至少为 1，games 至少为 2")
    weights = (
        BotWeights.load(arguments.weights)
        if arguments.weights.exists()
        else BotWeights()
    )
    train_layout(
        arguments.generations,
        arguments.games,
        arguments.seed,
        weights,
        arguments.output,
    )


if __name__ == "__main__":
    main()
