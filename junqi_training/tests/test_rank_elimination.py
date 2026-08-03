"""Tests for `OpponentKnowledge.eliminate_dead_ranks`.

The deduction was written after a language model produced it unprompted on a
real position: both enemy bombs accounted for, therefore the piece that traded
with our general *was* the general, therefore no survivor is one. The engine
had never done that -- it reasons one square at a time and never counts.

Soundness is the whole risk. A belief that wrongly excludes a rank is worse
than a loose one: the search will price an attack against a piece it has
"proven" cannot beat it, and walk into it. So the load-bearing test is not that
the deduction fires, it is that it is never wrong on a real game.
"""

import unittest

from junqi.arena import make_opening
from junqi.bot import BotWeights, HeuristicBot
from junqi.game import Game
from junqi.knowledge import OpponentKnowledge
from junqi.types import PIECE_COUNTS, Owner, PieceKind

WEIGHTS = BotWeights()


def _play(seed: int, plies: int):
    """Yield (game, knowledge) after each ply, tracking HUMAN's view of BOT."""
    game = Game(board=make_opening(seed), turn=Owner(seed % 2))
    players = {
        owner: HeuristicBot(WEIGHTS, seed=seed * 2 + int(owner)) for owner in Owner
    }
    knowledge = OpponentKnowledge(Owner.HUMAN)
    seen = 0
    for _ in range(plies):
        if game.over:
            return
        game.apply(players[game.turn].choose_move(game))
        for event in game.observations(Owner.HUMAN, seen):
            knowledge.observe(event)
        seen = len(game.records)
        knowledge.forget_missing(
            {p for p, piece in game.board.items() if piece.owner == Owner.BOT}
        )
        yield game, knowledge


class SoundnessTests(unittest.TestCase):
    def test_an_extinct_rank_really_has_no_survivors(self) -> None:
        """Soundness on real games.

        Note this test is currently *vacuous by measurement*, not by accident:
        over 60 full games (3,482 plies) the rule fired zero times, with and
        without the commander-alive signal. See
        `test_it_is_inert_at_this_rule_set`, which pins that as the known state
        so a future change that makes it fire shows up as a failing test rather
        than a silent behaviour change.
        """
        checked = 0
        for seed in range(20):
            for game, knowledge in _play(seed, 90):
                extinct = knowledge.eliminate_dead_ranks()
                alive = {
                    piece.kind
                    for piece in game.board.values()
                    if piece.owner == Owner.BOT
                }
                for kind in extinct:
                    self.assertNotIn(
                        kind,
                        alive,
                        f"seed {seed}: declared {kind.name} extinct while one lives",
                    )
                checked += 1
        self.assertGreater(checked, 500)

    def test_it_is_inert_at_this_rule_set(self) -> None:
        """The measured negative, pinned so a change to it is visible.

        Extinction needs `PIECE_COUNTS[R]` casualties pinned to one rank in one
        game. Measured over 60 games: 723 enemy deaths, of which only 6 had a
        singleton candidate set (0.8%) -- never enough. 30% of deaths are
        two-element `{R, BOMB}` sets, so the way to make this pay is a capacity
        argument over those (at most two of them can be bombs), not more
        singleton sources. Adding the commander-alive signal, which is what
        unlocked the chain for a language model on one position, did not change
        the count.
        """
        fired = 0
        for seed in range(20):
            for _, knowledge in _play(seed, 90):
                fired += bool(knowledge.eliminate_dead_ranks())
        self.assertEqual(
            fired,
            0,
            "the rule now fires -- re-measure its strength effect before "
            "shipping it, and update this test",
        )

    def test_a_live_piece_never_loses_its_true_rank_from_its_belief(self) -> None:
        """The invariant the search depends on: belief is loose, never wrong."""
        for seed in range(20):
            for game, knowledge in _play(seed, 90):
                knowledge.eliminate_dead_ranks()
                for square, kinds in knowledge.possible.items():
                    piece = game.board.get(square)
                    if piece is None or piece.owner != Owner.BOT:
                        continue
                    self.assertIn(
                        piece.kind,
                        kinds,
                        f"seed {seed} {square}: excluded the true rank",
                    )

    def test_it_never_empties_a_belief_set(self) -> None:
        for seed in range(15):
            for _, knowledge in _play(seed, 80):
                knowledge.eliminate_dead_ranks()
                for square, kinds in knowledge.possible.items():
                    self.assertTrue(kinds, f"seed {seed} {square}: belief emptied")

    def test_it_is_idempotent(self) -> None:
        for seed in range(5):
            for _, knowledge in _play(seed, 60):
                first = knowledge.eliminate_dead_ranks()
                before = dict(knowledge.possible)
                second = knowledge.eliminate_dead_ranks()
                self.assertEqual(first, second)
                self.assertEqual(before, knowledge.possible)


class DeductionTests(unittest.TestCase):
    def test_a_mine_that_trades_proves_a_bomb_died(self) -> None:
        """Nothing but a bomb can trade with a mine -- a singleton, provably.

        This is the step that unlocks the chain: two of these exhaust the
        enemy's bombs, after which any mutual destruction against our rank R
        pins the casualty to R exactly.
        """
        from junqi.knowledge import _dead_attacker_candidates

        self.assertEqual(
            _dead_attacker_candidates(PieceKind.MINE, 0), frozenset({PieceKind.BOMB})
        )

    def test_trading_with_our_colonel_leaves_colonel_or_bomb(self) -> None:
        from junqi.knowledge import _dead_attacker_candidates

        self.assertEqual(
            _dead_attacker_candidates(PieceKind.COLONEL, 0),
            frozenset({PieceKind.COLONEL, PieceKind.BOMB}),
        )

    def test_an_engineer_win_leaves_the_flag_in_and_so_pins_nothing(self) -> None:
        """Regression: excluding the flag here declared mines extinct wrongly.

        An engineer beats exactly the flag and mines. The last record of a won
        game is a flag capture, so dropping the flag makes that record read as
        a proven mine kill -- three of them "proved" mines extinct on seed 14
        with one still standing.
        """
        from junqi.knowledge import _dead_defender_candidates

        self.assertEqual(
            _dead_defender_candidates(PieceKind.ENGINEER, 1),
            frozenset({PieceKind.FLAG, PieceKind.MINE}),
        )

    def test_a_dead_attacker_is_never_a_mine_or_a_flag(self) -> None:
        from junqi.knowledge import _dead_attacker_candidates

        for defender in PieceKind:
            for outcome in (-1, 0):
                candidates = _dead_attacker_candidates(defender, outcome)
                self.assertNotIn(PieceKind.MINE, candidates)
                self.assertNotIn(PieceKind.FLAG, candidates)

    def test_enough_pinned_deaths_make_a_rank_extinct(self) -> None:
        knowledge = OpponentKnowledge(Owner.HUMAN)
        knowledge.possible = {(0, 0): frozenset({PieceKind.BOMB, PieceKind.COLONEL})}
        for _ in range(PIECE_COUNTS[PieceKind.BOMB]):
            knowledge.dead_enemy.append(frozenset({PieceKind.BOMB}))
        self.assertEqual(
            knowledge.eliminate_dead_ranks(), frozenset({PieceKind.BOMB})
        )
        self.assertEqual(knowledge.possible[(0, 0)], frozenset({PieceKind.COLONEL}))

    def test_exhausting_bombs_then_pins_a_trade_and_cascades(self) -> None:
        """The full chain, which is what the model did by hand."""
        knowledge = OpponentKnowledge(Owner.HUMAN)
        knowledge.possible = {(0, 0): frozenset({PieceKind.GENERAL, PieceKind.COLONEL})}
        # Two mines traded -> two bombs pinned -> bombs extinct.
        for _ in range(PIECE_COUNTS[PieceKind.BOMB]):
            knowledge.dead_enemy.append(frozenset({PieceKind.BOMB}))
        # Our general traded with something: general or bomb. Bombs are gone,
        # so it was their general -- and they only have one.
        knowledge.dead_enemy.append(frozenset({PieceKind.GENERAL, PieceKind.BOMB}))
        extinct = knowledge.eliminate_dead_ranks()
        self.assertEqual(extinct, frozenset({PieceKind.BOMB, PieceKind.GENERAL}))
        self.assertEqual(knowledge.possible[(0, 0)], frozenset({PieceKind.COLONEL}))

    def test_one_pinned_death_of_a_two_copy_rank_proves_nothing(self) -> None:
        knowledge = OpponentKnowledge(Owner.HUMAN)
        knowledge.possible = {(0, 0): frozenset({PieceKind.BOMB, PieceKind.COLONEL})}
        knowledge.dead_enemy.append(frozenset({PieceKind.BOMB}))
        self.assertEqual(knowledge.eliminate_dead_ranks(), frozenset())
        self.assertEqual(
            knowledge.possible[(0, 0)], frozenset({PieceKind.BOMB, PieceKind.COLONEL})
        )


class GatingTests(unittest.TestCase):
    def test_it_ships_off_so_current_behaviour_is_unchanged(self) -> None:
        self.assertEqual(BotWeights().use_rank_elimination, 0.0)

    def test_the_weight_reaches_both_belief_carrying_opponents(self) -> None:
        # A code-level `if` would apply to both sides of `compare()` and cancel
        # out (D022), so this has to be a coefficient the harness can vary.
        from dataclasses import replace

        from junqi.opponents import SelectiveBot
        from junqi.search_bot import SearchBot

        on = replace(WEIGHTS, use_rank_elimination=1.0)
        game = Game(board=make_opening(3), turn=Owner.HUMAN)
        for agent in (
            SearchBot(on, seed=1, samples=1, beam_width=2, reply_width=1),
            SelectiveBot(on, seed=1),
        ):
            move = agent.choose_move(game, Owner.HUMAN)
            self.assertIn(move, game.legal_moves(Owner.HUMAN))


if __name__ == "__main__":
    unittest.main()
