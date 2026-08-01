import random
import tempfile
import unittest
from collections import Counter
from dataclasses import asdict, replace
from pathlib import Path

from junqi.arena import compare, make_opening
from junqi.board import CAMPS, HEADQUARTERS
from junqi.bot import BotWeights, HeuristicBot
from junqi.cli import render
from junqi.opponents import AgentSpec, Pool, standard_pool
from junqi.web_export import parse as parse_web_weights
from junqi.deployment import (
    front_row,
    headquarters,
    load_deployment,
    random_deployment,
    save_deployment,
    strategic_deployment,
    swap_pieces,
    validate_deployment,
)
from junqi.game import Game, ObservedMove, battle_outcome
from junqi.knowledge import OpponentKnowledge
from junqi.search_bot import SearchBot
from junqi.types import Move, Owner, PIECE_COUNTS, Piece, PieceKind


class GameTests(unittest.TestCase):
    def test_new_game_has_complete_armies_and_empty_camps(self) -> None:
        game = Game.new(seed=7)
        for owner in Owner:
            pieces = [piece for piece in game.board.values() if piece.owner == owner]
            self.assertEqual(len(pieces), sum(PIECE_COUNTS.values()))
            self.assertEqual(len(pieces), 25)
            self.assertEqual(
                sum(piece.kind == PieceKind.FLAG for piece in pieces), 1
            )
        self.assertFalse(set(game.board) & CAMPS)

    def test_setup_is_reproducible(self) -> None:
        self.assertEqual(Game.new(seed=123).board, Game.new(seed=123).board)

    def test_battle_rules(self) -> None:
        self.assertEqual(battle_outcome(PieceKind.ENGINEER, PieceKind.MINE), 1)
        self.assertEqual(battle_outcome(PieceKind.COMMANDER, PieceKind.MINE), -1)
        self.assertEqual(battle_outcome(PieceKind.BOMB, PieceKind.COMMANDER), 0)
        self.assertEqual(battle_outcome(PieceKind.GENERAL, PieceKind.BRIGADIER), 1)
        self.assertEqual(battle_outcome(PieceKind.CAPTAIN, PieceKind.CAPTAIN), 0)

    def test_every_generated_move_can_be_applied(self) -> None:
        game = Game.new(seed=19)
        self.assertTrue(game.legal_moves())
        for move in game.legal_moves()[:20]:
            clone = game.clone()
            clone.apply(move)
            self.assertEqual(clone.move_count, 1)

    def test_turn_changes_after_quiet_move(self) -> None:
        game = Game.new(seed=4)
        quiet_move = next(
            move for move in game.legal_moves() if move.dst not in game.board
        )
        game.apply(quiet_move)
        self.assertEqual(game.turn, Owner.BOT)
        self.assertTrue(game.history)

    def test_headquarters_piece_cannot_move(self) -> None:
        game = Game.new(seed=9)
        for headquarters in HEADQUARTERS:
            if headquarters not in game.board:
                continue
            owner = game.board[headquarters].owner
            self.assertTrue(
                all(move.src != headquarters for move in game.legal_moves(owner))
            )

    def test_random_deployments_obey_position_rules(self) -> None:
        for seed in range(20):
            for owner in Owner:
                deployment = random_deployment(owner, random.Random(seed))
                self.assertEqual(validate_deployment(deployment, owner), [])
                bombs = [
                    position
                    for position, piece in deployment.items()
                    if piece.kind == PieceKind.BOMB
                ]
                self.assertTrue(
                    all(position[0] != front_row(owner) for position in bombs)
                )

    def test_invalid_swap_is_rejected_and_rolled_back(self) -> None:
        game = Game.new(seed=12)
        flag = game.flag_position(Owner.HUMAN)
        self.assertIsNotNone(flag)
        non_headquarters = next(
            position
            for position, piece in game.board.items()
            if piece.owner == Owner.HUMAN and position not in HEADQUARTERS
        )
        original = game.board.copy()
        with self.assertRaises(ValueError):
            swap_pieces(game.board, Owner.HUMAN, flag, non_headquarters)
        self.assertEqual(game.board, original)

    def test_battle_keeps_ranks_hidden(self) -> None:
        game = Game(
            board={
                (6, 0): Piece(Owner.HUMAN, PieceKind.GENERAL),
                (5, 0): Piece(Owner.BOT, PieceKind.BRIGADIER),
                (0, 1): Piece(Owner.BOT, PieceKind.FLAG),
            }
        )
        result = game.apply(Move((6, 0), (5, 0)))
        survivor = game.board[(5, 0)]
        self.assertFalse(survivor.revealed)
        self.assertNotIn("GENERAL", result)
        self.assertNotIn("BRIGADIER", result)
        self.assertIn("双方军衔保持隐藏", result)
        self.assertFalse(game.board[(0, 1)].revealed)

    def test_commander_death_reveals_own_flag(self) -> None:
        game = Game(
            board={
                (6, 0): Piece(Owner.HUMAN, PieceKind.BOMB),
                (5, 0): Piece(Owner.BOT, PieceKind.COMMANDER),
                (0, 1): Piece(Owner.BOT, PieceKind.FLAG),
                (11, 1): Piece(Owner.HUMAN, PieceKind.FLAG),
            }
        )
        result = game.apply(Move((6, 0), (5, 0)))
        self.assertTrue(game.board[(0, 1)].revealed)
        self.assertIn("军旗位置亮出：B1", result)
        self.assertIn("[旗]", render(game))

    def test_hidden_rank_does_not_change_heuristic_score(self) -> None:
        common = {
            (6, 0): Piece(Owner.HUMAN, PieceKind.GENERAL),
            (11, 1): Piece(Owner.HUMAN, PieceKind.FLAG),
            (0, 1): Piece(Owner.BOT, PieceKind.FLAG),
        }
        weak = Game({**common, (5, 0): Piece(Owner.BOT, PieceKind.LIEUTENANT)})
        strong = Game({**common, (5, 0): Piece(Owner.BOT, PieceKind.COMMANDER)})
        move = Move((6, 0), (5, 0))
        weak_score = HeuristicBot(BotWeights(), 3)._score(weak, move, Owner.HUMAN)
        strong_score = HeuristicBot(BotWeights(), 3)._score(
            strong, move, Owner.HUMAN
        )
        self.assertEqual(weak_score, strong_score)

    def test_search_bot_returns_legal_move(self) -> None:
        game = Game.new(seed=21, first=Owner.BOT)
        move = SearchBot(BotWeights(), seed=21, samples=1, beam_width=2).choose_move(
            game
        )
        self.assertIn(move, game.legal_moves())

    def test_bot_privately_infers_survivor_rank_range(self) -> None:
        game = Game(
            board={
                (6, 0): Piece(Owner.HUMAN, PieceKind.GENERAL),
                (5, 0): Piece(Owner.BOT, PieceKind.BRIGADIER),
                (0, 1): Piece(Owner.BOT, PieceKind.FLAG),
                (11, 1): Piece(Owner.HUMAN, PieceKind.FLAG),
            }
        )
        game.apply(Move((6, 0), (5, 0)))
        observation = game.observations(Owner.BOT)[0]
        knowledge = OpponentKnowledge(Owner.BOT)
        knowledge.observe(observation)
        self.assertEqual(
            knowledge.possible[(5, 0)],
            frozenset(
                {
                    PieceKind.COMMANDER,
                    PieceKind.GENERAL,
                    PieceKind.MAJOR_GENERAL,
                }
            ),
        )
        knowledge.observe(
            ObservedMove(
                Move((5, 0), (4, 0)),
                Owner.HUMAN,
                False,
                None,
                None,
            )
        )
        self.assertNotIn((5, 0), knowledge.possible)
        self.assertIn((4, 0), knowledge.possible)

    def test_layout_model_round_trip(self) -> None:
        layout = random_deployment(Owner.BOT, random.Random(44))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "layout.json"
            save_deployment(layout, path)
            loaded = load_deployment(path)
        self.assertEqual(layout, loaded)
        game = Game.new(seed=3, bot_deployment=loaded)
        self.assertEqual(
            {
                position: piece
                for position, piece in game.board.items()
                if piece.owner == Owner.BOT
            },
            loaded,
        )

    def test_flag_candidates_narrow_as_headquarters_fall(self) -> None:
        game = Game.new(seed=17)
        self.assertEqual(len(game.flag_candidates(Owner.BOT)), 2)
        decoy = next(
            position
            for position in game.flag_candidates(Owner.BOT)
            if game.board[position].kind != PieceKind.FLAG
        )
        game.board.pop(decoy)
        self.assertEqual(
            game.flag_candidates(Owner.BOT), [game.flag_position(Owner.BOT)]
        )

    def test_bot_storms_the_last_headquarters(self) -> None:
        game = Game(
            board={
                (0, 1): Piece(Owner.BOT, PieceKind.FLAG),
                (1, 1): Piece(Owner.HUMAN, PieceKind.CAPTAIN),
                (1, 2): Piece(Owner.HUMAN, PieceKind.MAJOR),
                (11, 1): Piece(Owner.HUMAN, PieceKind.FLAG),
            }
        )
        bot = HeuristicBot(BotWeights(noise=0.0), seed=1)
        strike = bot._score(game, Move((1, 1), (0, 1)), Owner.HUMAN)
        quiet = bot._score(game, Move((1, 2), (1, 3)), Owner.HUMAN)
        self.assertGreater(strike, quiet + 50)
        self.assertEqual(bot.choose_move(game, Owner.HUMAN), Move((1, 1), (0, 1)))

    def test_engineers_are_spent_on_mines_not_on_midboard_unknowns(self) -> None:
        weights = BotWeights(noise=0.0)

        def score(attacker: PieceKind, source: tuple[int, int], target: tuple[int, int]) -> float:
            game = Game(
                board={
                    source: Piece(Owner.HUMAN, attacker),
                    target: Piece(Owner.BOT, PieceKind.MAJOR),
                    (11, 1): Piece(Owner.HUMAN, PieceKind.FLAG),
                }
            )
            return HeuristicBot(weights, seed=2)._score(
                game, Move(source, target), Owner.HUMAN
            )

        # Rear rows are the only place a mine can be, and only engineers survive one.
        self.assertGreater(
            score(PieceKind.ENGINEER, (2, 0), (1, 0)),
            score(PieceKind.COLONEL, (2, 0), (1, 0)),
        )
        # Anywhere else an engineer loses to every rank in the game.
        self.assertLess(
            score(PieceKind.ENGINEER, (5, 0), (4, 0)),
            score(PieceKind.COLONEL, (5, 0), (4, 0)),
        )

    def test_strategic_deployments_are_legal_and_unpredictable(self) -> None:
        signatures = set()
        flag_squares = set()
        for seed in range(30):
            for owner in Owner:
                layout = strategic_deployment(owner, random.Random(seed))
                self.assertEqual(validate_deployment(layout, owner), [])
                occupants = {
                    position: layout[position].kind
                    for position in headquarters(owner)
                }
                decoys = [
                    kind
                    for kind in occupants.values()
                    if kind != PieceKind.FLAG
                ]
                self.assertEqual(len(decoys), 1)
                # A headquarters piece can never move, so it must be cheap.
                self.assertIn(
                    decoys[0], {PieceKind.LIEUTENANT, PieceKind.CAPTAIN}
                )
                if owner == Owner.BOT:
                    signatures.add(
                        tuple(
                            sorted(
                                (position, piece.kind.name)
                                for position, piece in layout.items()
                            )
                        )
                    )
                    flag_squares.add(
                        next(
                            position
                            for position, kind in occupants.items()
                            if kind == PieceKind.FLAG
                        )
                    )
        self.assertEqual(len(signatures), 30)
        self.assertEqual(len(flag_squares), 2)

    def test_determinize_does_not_copy_true_survivor_ranks(self) -> None:
        game = Game.new(seed=31)
        victims = [
            position
            for position, piece in game.board.items()
            if piece.owner == Owner.HUMAN and piece.kind == PieceKind.CAPTAIN
        ][:2]
        for position in victims:
            game.board.pop(position)
        bot = SearchBot(BotWeights(), seed=31)
        alive = sum(
            1 for piece in game.board.values() if piece.owner == Owner.HUMAN
        )
        multisets = set()
        for seed in range(12):
            world = bot._determinize(game, Owner.HUMAN, random.Random(seed))
            counts = Counter(
                piece.kind
                for piece in world.board.values()
                if piece.owner == Owner.HUMAN
            )
            self.assertEqual(sum(counts.values()), alive)
            self.assertEqual(counts[PieceKind.FLAG], 1)
            multisets.add(tuple(sorted(counts.items())))
        # Battles are anonymous: the bot must not know which ranks died.
        self.assertGreater(len(multisets), 1)

    def test_web_weights_are_in_sync_with_the_trained_model(self) -> None:
        root = Path(__file__).resolve().parents[1]
        generated = root / "web" / "lib" / "weights.ts"
        if not generated.exists():
            # The training subset is shipped to the compute box without `web/`
            # (see scripts/sync_remote.py); there is nothing to be in sync with.
            self.skipTest("no web engine in this checkout")
        model = BotWeights.load(root / "models" / "bot_weights.json")
        exported = parse_web_weights(generated.read_text(encoding="utf-8"))
        expected = asdict(model)
        self.assertEqual(
            set(exported),
            set(expected),
            "run `python -m junqi.web_export` after changing BotWeights",
        )
        for name, value in expected.items():
            self.assertAlmostEqual(exported[name], value, places=6, msg=name)

    def test_openings_depend_only_on_the_seed(self) -> None:
        self.assertEqual(make_opening(7), make_opening(7))
        self.assertNotEqual(make_opening(7), make_opening(8))
        for owner in Owner:
            self.assertEqual(validate_deployment(make_opening(7), owner), [])

    def test_sampled_worlds_do_not_depend_on_the_weights(self) -> None:
        """Common random numbers: the whole comparison design rests on this."""

        class Recording(SearchBot):
            def __init__(self, *args, **kwargs) -> None:
                super().__init__(*args, **kwargs)
                self.seen: list[int] = []

            def _rollout(self, game, move, owner, sample_seed):
                self.seen.append(sample_seed)
                return super()._rollout(game, move, owner, sample_seed)

        game = Game(board=make_opening(3), turn=Owner.BOT)
        plain = Recording(BotWeights(), seed=11, samples=4, beam_width=3)
        altered = Recording(
            replace(BotWeights(), forward=9.9, capture=-3.0, hq_strike=99.0),
            seed=11,
            samples=4,
            beam_width=3,
        )
        plain.choose_move(game, Owner.BOT)
        altered.choose_move(game, Owner.BOT)
        self.assertEqual(sorted(set(plain.seen)), sorted(set(altered.seen)))
        self.assertEqual(len(set(plain.seen)), 4)

    def test_every_pool_opponent_plays_legally(self) -> None:
        game = Game(board=make_opening(5), turn=Owner.HUMAN)
        pool = standard_pool()
        self.assertGreaterEqual(len(pool), 10)
        for spec in pool.specs:
            agent = spec.build(BotWeights(), seed=5)
            move = agent.choose_move(game, Owner.HUMAN)
            self.assertIn(move, game.legal_moves(Owner.HUMAN), spec.name)

    def test_an_anchored_pool_does_not_track_the_model_under_test(self) -> None:
        """Opponents must be identical on both sides of a comparison.

        Without an anchor, `play_match` builds a weight-driven opponent from
        the *subject's* weights, so a candidate is judged against a distorted
        copy of itself while the baseline is judged against a copy of itself --
        and the pool stops being a yardstick.
        """
        loose = {spec.name for spec in standard_pool().specs if spec.weights_path}
        self.assertEqual(loose, set())

        anchored = standard_pool(anchor="models/defaults.json")
        weight_driven = [
            spec for spec in anchored.specs if spec.kind in {"heuristic", "search"}
        ]
        self.assertTrue(weight_driven)
        for spec in weight_driven:
            self.assertEqual(spec.weights_path, "models/defaults.json", spec.name)
        # random / hqrush ignore weights entirely, so they need no anchor.
        for spec in anchored.specs:
            if spec.kind in {"random", "hqrush"}:
                self.assertIsNone(spec.weights_path, spec.name)

    def test_identical_weights_produce_a_zero_paired_difference(self) -> None:
        pool = Pool([AgentSpec("quiet", "heuristic", noise=0.0)])
        verdict = compare(BotWeights(), BotWeights(), pool, [1, 2], workers=1)
        self.assertEqual(verdict.candidate.games, verdict.incumbent.games)
        self.assertAlmostEqual(verdict.mean_difference, 0.0, places=12)
        self.assertFalse(verdict.significant)

    def test_only_a_c_e_files_cross_the_river(self) -> None:
        game = Game(
            board={
                (5, 0): Piece(Owner.HUMAN, PieceKind.LIEUTENANT),
                (5, 1): Piece(Owner.HUMAN, PieceKind.CAPTAIN),
                (5, 2): Piece(Owner.HUMAN, PieceKind.ENGINEER),
                (5, 3): Piece(Owner.HUMAN, PieceKind.MAJOR),
                (5, 4): Piece(Owner.HUMAN, PieceKind.COLONEL),
            }
        )
        moves = set(game.legal_moves(Owner.HUMAN))
        self.assertIn(Move((5, 0), (6, 0)), moves)
        self.assertNotIn(Move((5, 1), (6, 1)), moves)
        self.assertIn(Move((5, 2), (6, 2)), moves)
        self.assertNotIn(Move((5, 3), (6, 3)), moves)
        self.assertIn(Move((5, 4), (6, 4)), moves)


if __name__ == "__main__":
    unittest.main()
