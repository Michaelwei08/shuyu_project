"""Drive the real `play()` to a finish and check the replay it writes.

This test exists because of a specific failure. `--replay` was added, smoke
tested by calling `format_replay` directly, and shipped -- while `play()` itself
raised `NameError: name 'opening' is not defined` at the write, *after* the game
was over. A human played a full game against `oracle-perfect`, won it, and got a
traceback instead of a record. The formatter was tested; the path was not.

So the assertion here is deliberately end-to-end: patch `input`, let `play()`
run a whole game and write the file, then parse what landed on disk. Nothing is
called directly except `play`.
"""

from __future__ import annotations

import builtins
import contextlib
import io
import re
import tempfile
import unittest
from pathlib import Path

from junqi.cli import SYMBOLS, play
from junqi.game import Game
from junqi.opponents import HQRushBot
from junqi.types import Owner, PieceKind

#: The human is played by the opponent that charges the flag, so the game
#: actually ends. A random mover can wander for thousands of plies -- and
#: `play()` has no move cap, because a person gets bored and quits instead.
MAX_INPUTS = 600


def _run(
    opponent: str,
    seed: int,
    destination: Path,
    arrange: list[str] | None = None,
) -> tuple[str, dict[str, Game]]:
    """Play one whole game through `play()` and return the replay it wrote.

    ``arrange`` drives the manual deployment phase instead of `--auto-deploy`.
    That path matters on its own: the snapshot is deliberately taken *after*
    arranging, so a replay must disclose the board actually played rather than
    the pre-swap one, and only this branch can catch it being taken too early.
    """
    holder: dict[str, Game] = {}
    original = Game.new
    calls = {"count": 0}
    setup = list(arrange or [])

    def capture(*args, **kwargs):
        game = original(*args, **kwargs)
        holder["game"] = game
        # The layout before `arrange_player` touches it. The final board cannot
        # stand in for it -- pieces die, and the loser's flag is removed
        # outright -- so the pre-arrange copy is the only thing a disclosed
        # opening can be checked against.
        holder["pre"] = dict(game.board)
        return game

    def feed(_prompt: str = "") -> str:
        if setup:
            return setup.pop(0)
        calls["count"] += 1
        if calls["count"] > MAX_INPUTS:
            raise AssertionError(f"{opponent} game did not finish in {MAX_INPUTS} plies")
        game = holder["game"]
        # A fresh bot each turn: it keeps no state that matters here, and this
        # keeps the test independent of how HQRushBot is constructed.
        return str(HQRushBot(seed=calls["count"]).choose_move(game, Owner.HUMAN))

    Game.new = classmethod(lambda cls, **kwargs: capture(**kwargs))  # type: ignore[assignment]
    real_input = builtins.input
    builtins.input = feed
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            play(
                seed=seed,
                model_path=Path("models/bot_weights.json"),
                auto_deploy=arrange is None,
                opponent=opponent,
                replay_path=destination,
            )
    finally:
        builtins.input = real_input
        Game.new = original  # type: ignore[assignment]
    return destination.read_text(encoding="utf-8"), holder


class ReplayTests(unittest.TestCase):
    def test_a_full_game_writes_a_replay_that_parses(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            destination = Path(folder) / "game.txt"
            text, _ = _run("oracle-perfect", 20260803, destination)

        lines = text.split("\n")
        self.assertEqual(lines[0], "JQ/60 replay")
        self.assertRegex(lines[1], r"^result:     (human|bot) wins at ply \d+$")
        self.assertEqual(lines[2], "difficulty: oracle-perfect")
        self.assertEqual(lines[3], "engine:     python cli")
        self.assertRegex(lines[4], r"^weights:    [0-9a-f]{8}$")

        # The bug that made this test necessary: the opening snapshot.
        self.assertIn("bot opening (rows 1-6, disclosed now the game is over):", text)
        self.assertIn("your opening (rows 7-12):", text)

        moves = re.findall(
            r"^\s*(\d+) ([YB])  ([A-E]\d{1,2})-([A-E]\d{1,2})  \S+ → \S+(.*)$",
            text,
            re.MULTILINE,
        )
        self.assertGreater(len(moves), 10)
        # Plies are numbered from 1 with no gaps, and the sides alternate.
        for index, (number, side, _src, _dst, _note) in enumerate(moves, start=1):
            self.assertEqual(int(number), index)
        sides = {side for _, side, _, _, _ in moves}
        self.assertEqual(sides, {"Y", "B"})

        # The result line has to agree with the last move.
        ply_count = int(re.search(r"wins at ply (\d+)", lines[1]).group(1))
        self.assertEqual(ply_count, len(moves))

    def test_a_hand_deployed_game_discloses_the_board_actually_played(self) -> None:
        """The path a person uses, and the one the snapshot ordering is for."""
        with tempfile.TemporaryDirectory() as folder:
            destination = Path(folder) / "game.txt"
            # Row 10 is neither the front row nor a rear row, so it can hold
            # neither a mine nor the flag and any two of its squares can always
            # trade -- which keeps this test from depending on the layout the
            # seed happened to produce.
            text, holder = _run(
                "oracle-perfect",
                777,
                destination,
                arrange=["swap A10 C10", "done"],
            )
        rows = re.search(
            r"your opening \(rows 7-12\):\n     A  B  C  D  E\n((?:.*\n){6})", text
        )
        self.assertIsNotNone(rows, "the human opening block is missing")
        disclosed: dict[tuple[int, int], str] = {}
        for offset, line in enumerate(rows.group(1).rstrip("\n").split("\n")):
            cells = line[3:]
            for column in range(5):
                symbol = cells[column * 2 : column * 2 + 2].strip()
                if symbol and symbol != "·":
                    disclosed[(offset + 6, column)] = symbol

        left, right = (9, 0), (9, 2)
        expected = {
            square: SYMBOLS[piece.kind]
            for square, piece in holder["pre"].items()
            if piece.owner == Owner.HUMAN
        }
        expected[left], expected[right] = expected[right], expected[left]
        self.assertEqual(
            disclosed,
            expected,
            "the disclosed opening is the pre-swap board, so the snapshot is "
            "being taken before arrange_player instead of after it",
        )

    def test_both_openings_are_disclosed_in_full(self) -> None:
        """A replay whose openings are blank is the failure that started this."""
        with tempfile.TemporaryDirectory() as folder:
            destination = Path(folder) / "game.txt"
            text, _ = _run("search", 4242, destination)

        pieces = set(SYMBOLS.values())
        grids = re.findall(r"^\s*\d+((?: (?:\S|\s))+)$", text, re.MULTILINE)
        # 12 rows, each holding 5 cells; every non-camp cell is a real symbol.
        rows = [row for row in grids if len(row) == 10]
        self.assertEqual(len(rows), 12, "expected 12 disclosed board rows")
        occupied = sum(
            1
            for row in rows
            for cell in (row[index : index + 2] for index in range(0, 10, 2))
            if cell.strip() in pieces
        )
        # 25 pieces a side at the start, minus whatever died before the end.
        self.assertGreater(occupied, 20)
        self.assertIn(SYMBOLS[PieceKind.FLAG], text)


if __name__ == "__main__":
    unittest.main()
