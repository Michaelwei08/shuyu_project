"""Publish trained weights to the TypeScript engine.

`web/lib/bot.ts` used to carry its coefficients as inline literals, so training
never reached the browser game at all. It now imports `web/lib/weights.ts`,
which this module generates from `models/bot_weights.json`.

Run after every accepted model:

    python -m junqi.web_export

`tests/test_game.py::test_web_weights_are_in_sync` fails if you forget.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from dataclasses import asdict, fields
from pathlib import Path

from .bot import BotWeights

DEFAULT_MODEL = Path("models/bot_weights.json")
DEFAULT_TARGET = Path("web/lib/weights.ts")

HEADER = """\
// GENERATED FILE -- do not edit by hand.
// Written by `python -m junqi.web_export` from models/bot_weights.json, so the
// browser bot and the Python bot are driven by one trained model. Editing a
// number here will be overwritten and will fail the Python sync test.

export type BotWeights = {
"""


def fingerprint(weights: BotWeights) -> str:
    """Short digest so a pasted game replay says which model produced it."""
    payload = json.dumps(asdict(weights), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:8]


def render(weights: BotWeights) -> str:
    lines = [HEADER]
    for descriptor in fields(weights):
        lines.append(f"  {descriptor.name}: number;\n")
    lines.append("};\n\nexport const WEIGHTS: BotWeights = {\n")
    payload = asdict(weights)
    for descriptor in fields(weights):
        value = payload[descriptor.name]
        lines.append(f"  {descriptor.name}: {_number(value)},\n")
    lines.append("};\n\n")
    lines.append(f'export const WEIGHTS_FINGERPRINT = "{fingerprint(weights)}";\n')
    return "".join(lines)


def _number(value: float) -> str:
    if value == int(value):
        return str(int(value))
    return repr(round(value, 10))


def parse(source: str) -> dict[str, float]:
    """Read the values back out of the generated file (used by the sync test)."""
    body = source.split("export const WEIGHTS: BotWeights = {", 1)
    if len(body) != 2:
        raise ValueError("weights.ts is missing its WEIGHTS literal")
    entries = re.findall(r"(\w+):\s*(-?[\d.eE+-]+),", body[1])
    return {name: float(value) for name, value in entries}


VALUE_HEADER = """\
// GENERATED FILE -- do not edit by hand.
// Written by `python -m junqi.web_export` from models/value.json.
//
// REFERENCE_* below is a fixed position and the prediction Python produces for
// it. `lib/value.ts` reimplements the feature extractor, and the two must agree
// exactly or the exported coefficients mean different things in the two
// engines; `tests/bot-runtime.test.mjs` recomputes the reference and fails on
// any drift.
"""


def render_value(model: "ValueModel | None") -> str:
    """Emit the learned leaf evaluation plus a cross-engine drift check."""
    from .arena import make_opening
    from .game import Game
    from .types import Owner, format_position

    if model is None:
        return (
            VALUE_HEADER
            + "\nexport const VALUE_WEIGHTS: number[] = [];\n"
            + "export const REFERENCE_BOARD: Record<string, string> = {};\n"
            + "export const REFERENCE_MOVE_COUNT = 0;\n"
            + "export const REFERENCE_PREDICTION = 0.5;\n"
        )

    # Play into the position rather than scoring the opening. A starting board
    # has equal material, equal piece counts and symmetric reach, so most
    # features are 0 or identical for both sides and a drift check against it
    # passes even when the extractor is wrong -- which it duly did.
    from .opponents import HQRushBot
    from .search_bot import SearchBot

    game = Game(board=make_opening(20_260_801), turn=Owner.BOT)
    players = {
        Owner.BOT: SearchBot(BotWeights(), seed=5, samples=1, beam_width=3),
        Owner.HUMAN: HQRushBot(seed=6),
    }
    while game.move_count < 34 and not game.over:
        game.apply(players[game.turn].choose_move(game))
    prediction = model.predict(game, Owner.BOT)
    squares = {
        format_position(position): f"{piece.owner.name}:{piece.kind.name}"
        for position, piece in sorted(game.board.items())
    }
    lines = [VALUE_HEADER, "\nexport const VALUE_WEIGHTS: number[] = [\n"]
    for index in range(0, len(model.weights), 6):
        row = ", ".join(_number(v) for v in model.weights[index : index + 6])
        lines.append(f"  {row},\n")
    lines.append("];\n\nexport const REFERENCE_BOARD: Record<string, string> = {\n")
    for coordinate, descriptor in squares.items():
        lines.append(f'  "{coordinate}": "{descriptor}",\n')
    lines.append("};\n\n")
    lines.append(f"export const REFERENCE_MOVE_COUNT = {game.move_count};\n")
    lines.append(f"export const REFERENCE_PREDICTION = {prediction!r};\n")
    return "".join(lines)


def export(model: Path, target: Path) -> BotWeights:
    weights = BotWeights.load(model)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(render(weights), encoding="utf-8", newline="\n")

    from .value import DEFAULT_PATH, ValueModel

    value_model = ValueModel.load(DEFAULT_PATH) if DEFAULT_PATH.exists() else None
    (target.parent / "value-weights.ts").write_text(
        render_value(value_model), encoding="utf-8", newline="\n"
    )
    return weights


def main() -> None:
    parser = argparse.ArgumentParser(description="把训练好的权重导出到网页引擎")
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--target", type=Path, default=DEFAULT_TARGET)
    arguments = parser.parse_args()
    if not arguments.model.exists():
        raise SystemExit(f"找不到模型文件：{arguments.model}")
    weights = export(arguments.model, arguments.target)
    print(
        f"已导出 {len(fields(weights))} 项权重 -> {arguments.target.resolve()}"
    )


if __name__ == "__main__":
    main()
