"""A piece that ends on a headquarters is frozen; the score must know it.

Immobility attaches to the *square*, not to the piece deployed there, so
winning a headquarters probe entombs the winner. Before `eval_hq_entomb` the
headquarters branch was rank-independent -- `hq_strike` alone -- so the
heuristic could not tell an engineer probe from a commander probe, while the
rollout preferred the commander because it survives. Measured over 320 pool
games, COMMANDER was entombed 41 times against ENGINEER's 1.
"""

from __future__ import annotations

import unittest

from junqi.board import HEADQUARTERS
from junqi.bot import BotWeights, HeuristicBot, _piece_value
from junqi.game import Game
from junqi.types import Move, Owner, Piece, PieceKind


def _probe_position() -> tuple[Game, tuple[int, int], tuple[int, int]]:
    """A board where the human may attack either of the bot's headquarters.

    Built by hand rather than from a deployment, because what matters is only
    that both headquarters are occupied (so neither is a certain flag) and that
    the attacker stands next to one.
    """
    board: dict[tuple[int, int], Piece] = {
        (0, 1): Piece(Owner.BOT, PieceKind.FLAG),
        (0, 3): Piece(Owner.BOT, PieceKind.CAPTAIN),
        # The bot needs a flag *and* a second headquarters occupant for the
        # candidate list to hold two squares.
        (11, 1): Piece(Owner.HUMAN, PieceKind.FLAG),
    }
    approach = (1, 3)
    return Game(board=board, turn=Owner.HUMAN), approach, (0, 3)


class EntombTests(unittest.TestCase):
    def test_the_square_freezes_whatever_lands_on_it(self) -> None:
        game, approach, headquarters = _probe_position()
        self.assertIn(headquarters, HEADQUARTERS)
        game.board[headquarters] = Piece(Owner.HUMAN, PieceKind.GENERAL)
        moves = [
            move
            for move in game.legal_moves(Owner.HUMAN)
            if move.src == headquarters
        ]
        self.assertEqual(moves, [], "a piece on a headquarters must be frozen")

    def test_shipped_weights_leave_the_term_off(self) -> None:
        """It ships at 0 until it has a paired p-value."""
        self.assertEqual(BotWeights().eval_hq_entomb, 0.0)

    def test_off_by_default_the_branch_cannot_tell_the_ranks_apart(self) -> None:
        """The defect the term exists to fix, pinned so it stays fixed."""
        game, approach, headquarters = _probe_position()
        bot = HeuristicBot(BotWeights(), seed=1)
        bot.weights.noise = 0.0
        scores = {}
        for kind in (PieceKind.COMMANDER, PieceKind.ENGINEER):
            game.board[approach] = Piece(Owner.HUMAN, kind)
            scores[kind] = bot._score(
                game, Move(approach, headquarters), Owner.HUMAN, quick=True
            )
        self.assertAlmostEqual(scores[PieceKind.COMMANDER], scores[PieceKind.ENGINEER])

    def test_on_it_prefers_the_cheap_prober(self) -> None:
        game, approach, headquarters = _probe_position()
        weights = BotWeights()
        weights.noise = 0.0
        weights.eval_hq_entomb = 1.0
        bot = HeuristicBot(weights, seed=1)
        scores = {}
        for kind in (PieceKind.COMMANDER, PieceKind.ENGINEER):
            game.board[approach] = Piece(Owner.HUMAN, kind)
            scores[kind] = bot._score(
                game, Move(approach, headquarters), Owner.HUMAN, quick=True
            )
        self.assertLess(
            scores[PieceKind.COMMANDER],
            scores[PieceKind.ENGINEER],
            "the expensive prober must be the less attractive one",
        )
        gap = scores[PieceKind.ENGINEER] - scores[PieceKind.COMMANDER]
        expected = _piece_value(PieceKind.COMMANDER) - _piece_value(PieceKind.ENGINEER)
        self.assertAlmostEqual(gap, expected)

    def test_a_certain_flag_is_never_penalised(self) -> None:
        """Being frozen is free on the move that ends the game."""
        game, approach, headquarters = _probe_position()
        # Vacate the other headquarters, so the one being attacked is the only
        # candidate left and therefore certainly the flag -- the case the
        # penalty must not touch.
        del game.board[(0, 1)]
        game.board[headquarters] = Piece(Owner.BOT, PieceKind.FLAG)
        weights = BotWeights()
        weights.noise = 0.0
        plain = HeuristicBot(weights, seed=1)
        weights_on = BotWeights()
        weights_on.noise = 0.0
        weights_on.eval_hq_entomb = 1.0
        loud = HeuristicBot(weights_on, seed=1)
        game.board[approach] = Piece(Owner.HUMAN, PieceKind.COMMANDER)
        move = Move(approach, headquarters)
        self.assertEqual(game.flag_candidates(Owner.BOT), [headquarters])
        self.assertAlmostEqual(
            plain._score(game, move, Owner.HUMAN, quick=True),
            loud._score(game, move, Owner.HUMAN, quick=True),
        )


if __name__ == "__main__":
    unittest.main()
