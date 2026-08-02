import json
import random
import tempfile
import unittest
from collections import Counter
from dataclasses import asdict, fields, replace
from pathlib import Path

from junqi.arena import compare, make_opening
from junqi.board import CAMPS, HEADQUARTERS
from junqi.bot import PRIOR_BATTLE, BotWeights, HeuristicBot
from junqi.cli import render
from junqi.opponents import AgentSpec, Pool, SelectiveBot, standard_pool
from junqi.web_export import parse as parse_web_weights
from junqi.deployment import (
    SCREEN_MINE_CAP,
    rear_rows,
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
from junqi.value import (
    BASE_WIDTH,
    WIDTH,
    ValueModel,
    features as value_features,
    load_default as load_value_model,
)


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

    def test_a_revealed_flag_proves_that_commander_is_dead(self) -> None:
        """The user's deduction: their commander is gone, so whatever beat our
        major general can only be the general."""
        bot = SearchBot(BotWeights(), seed=4)
        alive = bot._sample_survivors(24, Counter(), random.Random(1))
        self.assertIn(PieceKind.COMMANDER, alive)
        gone = bot._sample_survivors(
            24, Counter(), random.Random(1), commander_dead=True
        )
        self.assertNotIn(PieceKind.COMMANDER, gone)

        board = {
            (0, 1): Piece(Owner.BOT, PieceKind.FLAG, True),  # revealed by death
            (0, 3): Piece(Owner.BOT, PieceKind.CAPTAIN),
            (11, 1): Piece(Owner.HUMAN, PieceKind.FLAG),
        }
        game = Game(board=board, turn=Owner.HUMAN)
        self.assertTrue(bot._commander_dead(game, Owner.BOT))
        self.assertFalse(bot._commander_dead(game, Owner.HUMAN))

    def test_surviving_mines_are_not_killed_off_by_the_casualty_estimate(self) -> None:
        """Four rear pieces and no mine taken means three mines and the flag."""
        bot = SearchBot(BotWeights(), seed=9)
        # Squeeze the opponent down to four hidden pieces with nothing proven.
        for seed in range(20):
            alive = Counter(
                bot._sample_survivors(4, Counter(), random.Random(seed))
            )
            self.assertEqual(alive[PieceKind.MINE], 3, "a mine was invented dead")
            self.assertEqual(alive[PieceKind.FLAG], 1)

        # Once our engineer proves one dead, the estimate must let it go.
        proven = Counter({PieceKind.MINE: 1})
        alive = Counter(
            bot._sample_survivors(4, Counter(), random.Random(1), destroyed=proven)
        )
        self.assertEqual(alive[PieceKind.MINE], 2)

    def test_live_mines_cannot_outnumber_the_rear_squares_holding_them(self) -> None:
        """A mine never moves, so it must be standing on one of its own rear
        squares. Only an engineer win proves a mine dead, and a bomb trades with
        every rank so a bomb kill is unprovable -- which left the estimate
        carrying all three mines forever. Occupancy bounds it without reading a
        single hidden rank.
        """
        bot = SearchBot(BotWeights(), seed=9)
        for slots in range(4):
            alive = Counter(
                bot._sample_survivors(
                    8, Counter(), random.Random(slots), rear_slots=slots
                )
            )
            self.assertLessEqual(alive[PieceKind.MINE], slots, f"slots={slots}")
        # Unbounded is still the old behaviour, so the estimate stays generous
        # whenever the caller has nothing to say about occupancy.
        loose = Counter(bot._sample_survivors(8, Counter(), random.Random(0)))
        self.assertEqual(loose[PieceKind.MINE], 3)

    def test_determinize_never_puts_a_mine_outside_the_rear_rows(self) -> None:
        """The failure the occupancy cap exists to prevent.

        With more believed-live mines than legal rear squares,
        `_assign_constrained` cannot place them, exhausts its attempts and falls
        through to an unconstrained shuffle -- which used to scatter mines
        across midfield, where the rules say they can never be.
        """
        board = {
            (0, 1): Piece(Owner.BOT, PieceKind.FLAG),
            (11, 1): Piece(Owner.HUMAN, PieceKind.FLAG),
        }
        # One rear square left, but plenty of midfield pieces to place.
        board[(10, 0)] = Piece(Owner.HUMAN, PieceKind.MINE)
        for column, row in enumerate([7, 7, 7, 7]):
            board[(row, column)] = Piece(Owner.HUMAN, PieceKind.CAPTAIN)
        for column in range(4):
            board[(5, column)] = Piece(Owner.BOT, PieceKind.CAPTAIN)
        game = Game(board=board)

        bot = SearchBot(BotWeights(), seed=3)
        rear = rear_rows(Owner.HUMAN)
        for seed in range(25):
            world = bot._determinize(game, Owner.HUMAN, random.Random(seed))
            stray = [
                position
                for position, piece in world.board.items()
                if piece.owner == Owner.HUMAN
                and piece.kind == PieceKind.MINE
                and position[0] not in rear
            ]
            self.assertEqual(stray, [], f"seed {seed}: mine outside the rear rows")

    def test_the_modelled_opponent_can_be_given_beliefs_about_our_army(self) -> None:
        """`_rollout` built its reply bot and never assigned `.knowledge`, so
        every sampled world assumed an opponent with no deductions while the bot
        itself ran a full `OpponentKnowledge`. `reply_insight` ships at 0, which
        reproduces that exactly; the mechanism has to work so it can be measured.
        """
        game = Game(board=make_opening(5), turn=Owner.BOT)
        blind = SearchBot(BotWeights(), seed=5)
        self.assertIsNone(blind._reply_belief(game, Owner.BOT, 77))

        seeing = SearchBot(replace(BotWeights(), reply_insight=1.0), seed=5)
        belief = seeing._reply_belief(game, Owner.BOT, 77)
        assert belief is not None
        own = {
            position
            for position, piece in game.board.items()
            if piece.owner == Owner.BOT
        }
        self.assertEqual(set(belief), own)
        # It may only ever look at our own pieces -- never a hidden enemy rank.
        for position, kinds in belief.items():
            self.assertEqual(kinds, frozenset({game.board[position].kind}))
        self.assertEqual(
            seeing._reply_belief(game, Owner.BOT, 77), belief, "must be seed-stable"
        )

    def test_an_engineer_win_proves_a_mine_died(self) -> None:
        knowledge = OpponentKnowledge(Owner.BOT)
        win = ObservedMove(
            move=Move((1, 0), (0, 0)), attacker_owner=Owner.BOT, had_battle=True,
            outcome=1, own_kind=PieceKind.ENGINEER,
        )
        knowledge.observe(win)
        self.assertEqual(knowledge.destroyed[PieceKind.MINE], 1)
        # A colonel winning says only "something weaker" -- not nameable.
        knowledge.observe(replace(win, own_kind=PieceKind.COLONEL))
        self.assertEqual(knowledge.destroyed[PieceKind.MINE], 1)

    def test_a_blind_attack_is_priced_by_odds_not_by_attacker_value(self) -> None:
        """The mispricing that lost nine of ten games on 2026-07-31.

        `unknown_risk` scaled the penalty by what the attacker is *worth*, but
        what decides an attack on a square we know nothing about is how likely
        it is to *lose* -- and across the rank order those run in opposite
        directions. The old term therefore scored the commander lowest (+1.80)
        and the engineer highest (+2.53) for the very same blind attack, while
        the replayed games had the commander at 3W/1T/0L and the engineer at
        1W/2T/5L. The bot was being paid to probe with the pieces that lose.
        """
        weights = BotWeights(noise=0.0)

        def score(attacker: PieceKind) -> float:
            game = Game(
                board={
                    (5, 0): Piece(Owner.BOT, attacker),
                    # Row 7 is neither rear rows nor a headquarters, so this is
                    # the genuinely-unknown branch rather than mine country.
                    (6, 0): Piece(Owner.HUMAN, PieceKind.MAJOR),
                    (0, 1): Piece(Owner.BOT, PieceKind.FLAG),
                    (11, 1): Piece(Owner.HUMAN, PieceKind.FLAG),
                }
            )
            return HeuristicBot(weights, seed=7)._score(
                game, Move((5, 0), (6, 0)), Owner.BOT
            )

        # Engineers have their own branch (`engineer_waste`); everything else
        # runs through the prior, and must come out ordered by rank.
        ladder = [
            PieceKind.COMMANDER,
            PieceKind.GENERAL,
            PieceKind.MAJOR_GENERAL,
            PieceKind.BRIGADIER,
            PieceKind.COLONEL,
            PieceKind.MAJOR,
            PieceKind.CAPTAIN,
            PieceKind.LIEUTENANT,
        ]
        scores = [score(kind) for kind in ladder]
        for stronger, weaker, high, low in zip(
            ladder, ladder[1:], scores, scores[1:], strict=False
        ):
            self.assertGreater(high, low, f"{stronger.name} vs {weaker.name}")
        # The sign has to flip somewhere, or this is only a rescaling: a
        # lieutenant walking into an unknown must be worse than not doing it.
        self.assertLess(scores[-1], score(PieceKind.COMMANDER) - 10.0)

    def test_the_blind_prior_excludes_ranks_that_cannot_be_there(self) -> None:
        # A flag never leaves a headquarters and a mine never leaves the rear
        # two rows, so neither may lift an engineer's odds on a midboard square.
        self.assertLess(PRIOR_BATTLE[PieceKind.ENGINEER], 0.0)
        self.assertGreater(PRIOR_BATTLE[PieceKind.COMMANDER], 0.0)
        # A bomb trades with every rank, so its expectation is exactly zero.
        self.assertAlmostEqual(PRIOR_BATTLE[PieceKind.BOMB], 0.0)

    def test_defender_supply_is_charged_before_the_breach_lands(self) -> None:
        """The supply term measured dead (+0.0008 +/- 0.0097 over 806 games) and
        now ships at 0, but the mechanism has to keep working or the ablation
        that retired it could not be re-run. So this switches it on explicitly
        rather than relying on the default.
        """
        weights = replace(
            BotWeights(),
            eval_hq_supply=8.0,
            # Leave only material and the supply term standing, so the two
            # boards below differ by exactly one thing.
            eval_hq_guard=0.0,
            eval_mobility=0.0,
            eval_commander=0.0,
            eval_hq_attack=0.0,
            eval_hq_attack_certain=0.0,
            eval_hq_defense=0.0,
            eval_hq_defense_certain=0.0,
            eval_hq_breach=0.0,
            eval_immobilize=0.0,
        )
        bot = SearchBot(weights, seed=5)

        def value(defenders: tuple[tuple[int, int], ...]) -> float:
            board = {
                (0, 1): Piece(Owner.BOT, PieceKind.FLAG),
                (1, 3): Piece(Owner.HUMAN, PieceKind.CAPTAIN),  # two moves out
                (11, 1): Piece(Owner.HUMAN, PieceKind.FLAG),
            }
            for square in defenders:
                board[square] = Piece(Owner.BOT, PieceKind.CAPTAIN)
            return bot._state_value(Game(board=board), Owner.BOT)

        # Same two pieces, same material: held on the flag's approaches, or
        # sent away where they can never answer.
        home = value(((0, 0), (0, 2)))
        away = value(((8, 0), (8, 4)))
        self.assertGreater(home, away, "supply term is switched off")
        self.assertAlmostEqual(home - away, 2 * weights.eval_hq_supply, places=6)

    def test_the_flag_screen_cap_is_honoured(self) -> None:
        """The cap is a real knob, which is how the harness caught a bad idea.

        Capping the screen at 2 was measured at -0.0985 +/- 0.0217 over 806
        paired games, so the shipped policy is back to a full seal. This asserts
        the *plumbing*, not the policy: a lower cap must genuinely keep mines off
        the flag's doors, or `--screen-cap` would silently compare nothing.
        """
        def screen_mines(owner: Owner, seed: int, cap: int | None) -> int:
            layout = strategic_deployment(owner, random.Random(seed), screen_cap=cap)
            flag = next(
                position
                for position, piece in layout.items()
                if piece.kind == PieceKind.FLAG
            )
            return sum(
                1
                for row, column in ((1, 0), (-1, 0), (0, 1), (0, -1))
                if (square := (flag[0] + row, flag[1] + column)) in layout
                and layout[square].kind == PieceKind.MINE
            )

        for seed in range(40):
            for owner in Owner:
                for cap in (1, 2, 3):
                    self.assertLessEqual(
                        screen_mines(owner, seed, cap), cap, f"seed {seed} cap {cap}"
                    )
                self.assertLessEqual(
                    screen_mines(owner, seed, None), SCREEN_MINE_CAP, f"seed {seed}"
                )
        # A cap of 1 must actually leave doors open, or the knob is inert.
        capped = sum(screen_mines(Owner.BOT, seed, 1) for seed in range(40))
        sealed = sum(screen_mines(Owner.BOT, seed, 3) for seed in range(40))
        self.assertLess(capped, sealed)

    def test_the_screen_cap_changes_only_the_subject_army(self) -> None:
        """What makes the deployment change measurable by the paired harness.

        Each side draws from its own stream, so raising the bot's cap cannot
        shift how many draws are left for the human -- otherwise a
        capped-vs-sealed comparison would be measuring a different opponent as
        well as a different screen.
        """
        sealed = make_opening(11, {Owner.BOT: 3})
        capped = make_opening(11, {Owner.BOT: 2})

        def army(board, owner):
            return {
                position: piece
                for position, piece in board.items()
                if piece.owner == owner
            }

        self.assertEqual(army(sealed, Owner.HUMAN), army(capped, Owner.HUMAN))
        self.assertNotEqual(army(sealed, Owner.BOT), army(capped, Owner.BOT))
        for owner in Owner:
            self.assertEqual(validate_deployment(sealed, owner), [])
            self.assertEqual(validate_deployment(capped, owner), [])

    def test_the_selective_opponent_declines_the_fights_it_should_lose(self) -> None:
        """The policy class the pool was missing.

        Every weight-driven opponent shares `capture` and `blind_battle` with
        the candidate, and `material` doubles `capture`, so the whole pool
        over-attacks in the same direction and a mispriced attack term nets to
        zero in self-play. `hqrush` is the only structurally different opponent
        and it carries the opposite bias, refusing captures outright.
        """
        board = {
            (6, 0): Piece(Owner.HUMAN, PieceKind.LIEUTENANT),
            (5, 0): Piece(Owner.BOT, PieceKind.MAJOR),
            (0, 1): Piece(Owner.BOT, PieceKind.FLAG),
            (11, 1): Piece(Owner.HUMAN, PieceKind.FLAG),
        }
        game = Game(board=dict(board))
        agent = SelectiveBot(BotWeights(noise=0.0), seed=3)
        attack = Move((6, 0), (5, 0))
        chosen = agent.choose_move(game, Owner.HUMAN)
        self.assertIn(chosen, game.legal_moves(Owner.HUMAN))
        self.assertNotEqual(chosen, attack, "a lieutenant should not probe blind")

        # Same square, same ignorance, a rank that wins the exchange 86% of the
        # time: now it is worth taking.
        strong = Game(
            board={**board, (6, 0): Piece(Owner.HUMAN, PieceKind.COMMANDER)}
        )
        bold = SelectiveBot(BotWeights(noise=0.0), seed=3)
        bold._update_knowledge(strong, Owner.HUMAN)
        self.assertTrue(bold._acceptable(strong, attack, Owner.HUMAN, []))
        agent._update_knowledge(game, Owner.HUMAN)
        self.assertFalse(agent._acceptable(game, attack, Owner.HUMAN, []))

    def test_the_learned_value_term_is_inert_until_it_is_switched_on(self) -> None:
        """`eval_value_scale` ships at 0, so a fitted `models/value.json` sitting
        on disk must not change a single evaluation until a paired run says it
        should. That is what lets the whole learned evaluation be A/B'd as a
        weight instead of as a code switch."""
        game = Game(board=make_opening(3), turn=Owner.BOT)
        plain = SearchBot(BotWeights(), seed=1)
        self.assertEqual(BotWeights().eval_value_scale, 0.0)

        model = ValueModel(weights=[0.25] * WIDTH)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "value.json"
            model.save(path)
            loaded = ValueModel.load(path)
        self.assertEqual(loaded.weights, model.weights)

        off = plain._state_value(game, Owner.BOT)
        loud = SearchBot(
            replace(BotWeights(), eval_value_scale=50.0), seed=1
        )._state_value(game, Owner.BOT)
        if load_value_model() is None:
            self.skipTest("no models/value.json in this checkout")
        self.assertNotEqual(off, loud, "the term is wired but does nothing")

    def test_value_features_are_bounded_and_phase_separated(self) -> None:
        """Only the active phase's block is populated, which is what lets one
        linear model hold different opinions early and late -- composition
        measured worse than useless in the opening and dominant after ply 45."""
        game = Game(board=make_opening(8), turn=Owner.BOT)
        self.assertEqual(WIDTH, BASE_WIDTH * 3)

        seen = set()
        for move_count, expected in ((0, 0), (30, 1), (90, 2)):
            game.move_count = move_count
            vector = value_features(game, Owner.BOT)
            self.assertEqual(len(vector), WIDTH)
            active = {index // BASE_WIDTH for index, x in enumerate(vector) if x}
            self.assertEqual(active, {expected}, f"ply {move_count}")
            seen |= active
            # Everything is a ratio or a one-hot; an unbounded feature would let
            # one term swamp the dot product on an unusual board.
            self.assertTrue(all(-2.0 <= x <= 2.0 for x in vector))
        self.assertEqual(seen, {0, 1, 2})

    def test_value_model_rejects_a_stale_width(self) -> None:
        """Adding a feature invalidates every fitted model. Loading one anyway
        would silently misalign every coefficient, so it has to be an error."""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "value.json"
            ValueModel(weights=[0.0] * (WIDTH - 1)).save(path)
            with self.assertRaises(ValueError):
                ValueModel.load(path)

    def test_the_anchor_pins_every_coefficient(self) -> None:
        """The judging pool's yardstick must not move when a field is added.

        `BotWeights.load` deliberately tolerates missing keys so archived models
        still load -- but applied to the anchor that silently resolves every
        newly added coefficient to *today's* default, so the yardstick shifts on
        any commit that adds a field, with the file untouched. It had drifted
        into running `unknown_risk` and `blind_battle` simultaneously: both
        blind-attack penalties at once, which ten of thirteen pool opponents
        then used. Rerun `python scripts/rebuild_anchor.py --write` if this
        fails, and re-baseline afterwards.
        """
        anchor = Path(__file__).resolve().parents[1] / "models" / "defaults.json"
        if not anchor.exists():
            self.skipTest("no anchor in this checkout")
        stored = json.loads(anchor.read_text(encoding="utf-8"))
        expected = {descriptor.name for descriptor in fields(BotWeights)}
        self.assertEqual(
            set(stored),
            expected,
            "models/defaults.json is missing fields, so the pool silently "
            "tracks the dataclass -- run scripts/rebuild_anchor.py --write",
        )
        # The two blind-attack pricings are alternatives, never both at once.
        loaded = BotWeights.load(anchor)
        self.assertFalse(
            loaded.unknown_risk > 0 and loaded.blind_battle > 0,
            "the anchor prices a blind attack twice over",
        )

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
