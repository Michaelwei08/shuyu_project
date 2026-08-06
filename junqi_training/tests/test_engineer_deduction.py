"""Reading the opponent's engineer off a railway corner.

`engineer_expose` has always priced *our* side of this leak: only engineers turn
corners on the railway, so making one announces the piece. `OpponentKnowledge`
threw the other side away -- it learns from battle outcomes alone, so a quiet
move updated nothing at all. This is the one rank fact a quiet move carries, and
it is a certainty rather than a prior.

Measured over 160 pool games before the deduction existed: 1.67 such moves per
game, in 89% of games, and 267 of 267 were genuinely engineers.
"""

from __future__ import annotations

import random
import unittest

from junqi.arena import make_opening
from junqi.board import engineer_only_move, road_neighbors
from junqi.bot import BotWeights
from junqi.game import Game
from junqi.knowledge import MOVABLE_KINDS, OpponentKnowledge
from junqi.types import Owner, PieceKind


def _find_quiet_rail_turn(limit: int = 200):
    """Play until an engineer makes a corner turn onto an empty square."""
    for seed in range(700_000, 700_000 + limit):
        game = Game(board=make_opening(seed), turn=Owner.HUMAN)
        rng = random.Random(seed)
        for _ in range(150):
            if game.over:
                break
            moves = game.legal_moves()
            turns = [
                move
                for move in moves
                if game.board[move.src].kind == PieceKind.ENGINEER
                and move.dst not in game.board
                and engineer_only_move(move.src, move.dst, set(game.board))
            ]
            chosen = turns[0] if turns else rng.choice(moves)
            mover = game.turn
            game.apply(chosen)
            if turns:
                return game, mover, chosen
    raise AssertionError("no quiet rail turn found; the search is broken")


class PredicateTests(unittest.TestCase):
    def test_a_road_step_is_never_engineer_only(self) -> None:
        board = make_opening(5)
        occupied = set(board)
        for square in list(board)[:40]:
            for neighbour in road_neighbors(square):
                self.assertFalse(engineer_only_move(square, neighbour, occupied))

    def test_the_record_stamps_it_from_the_pre_move_board(self) -> None:
        """It cannot be recomputed later: the rails clear as pieces move."""
        game, mover, move = _find_quiet_rail_turn()
        record = game.records[-1]
        self.assertEqual(record.move, move)
        self.assertTrue(record.engineer_only)
        self.assertEqual(record.attacker.kind, PieceKind.ENGINEER)
        # And it reaches the observation stream both sides read.
        observed = game.observations(mover.other)[-1]
        self.assertTrue(observed.engineer_only)


class DeductionTests(unittest.TestCase):
    def test_shipped_weights_leave_it_off(self) -> None:
        self.assertEqual(BotWeights().use_engineer_deduction, 0.0)

    def test_off_the_belief_stays_wide_and_on_it_collapses(self) -> None:
        game, mover, move = _find_quiet_rail_turn()
        watcher = mover.other
        beliefs = {}
        for enabled in (False, True):
            knowledge = OpponentKnowledge(watcher)
            knowledge.deduce_engineers = enabled
            for event in game.observations(watcher):
                knowledge.observe(event)
            beliefs[enabled] = knowledge.possible.get(move.dst)

        self.assertEqual(beliefs[False], MOVABLE_KINDS)
        self.assertEqual(beliefs[True], frozenset({PieceKind.ENGINEER}))
        # The deduction is not merely narrow, it is right.
        self.assertEqual(game.board[move.dst].kind, PieceKind.ENGINEER)

    def test_a_battle_move_is_unaffected(self) -> None:
        """Battles already carried information; this must not disturb them."""
        game = Game(board=make_opening(11), turn=Owner.HUMAN)
        for _ in range(80):
            if game.over:
                break
            game.apply(random.Random(game.move_count).choice(game.legal_moves()))
        for watcher in Owner:
            plain, loud = OpponentKnowledge(watcher), OpponentKnowledge(watcher)
            loud.deduce_engineers = True
            for event in game.observations(watcher):
                plain.observe(event)
                loud.observe(event)
            for square, kinds in plain.possible.items():
                if loud.possible.get(square) == kinds:
                    continue
                # Any square that differs must be one the deduction pinned.
                self.assertEqual(
                    loud.possible[square], frozenset({PieceKind.ENGINEER})
                )


if __name__ == "__main__":
    unittest.main()
