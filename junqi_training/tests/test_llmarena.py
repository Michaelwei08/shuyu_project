"""Tests for the LLM arena scaffolding.

The one that matters is
``test_the_prompt_is_invariant_under_permuting_hidden_ranks``: it is the text
side of the engine's central invariant, and it is the reason this package can
be trusted to build prompts at all.
"""

import json
import pickle
import random
import tempfile
import unittest
from pathlib import Path

from junqi.arena import Job, make_opening, play_match
from junqi.bot import BotWeights, HeuristicBot
from junqi.game import Game
from junqi.opponents import AgentSpec
from junqi.search_bot import SearchBot
from junqi.types import Owner, Piece, PieceKind, format_position

from llmarena.agent import (
    LLMAgent,
    first_illegal_ply,
    parse_response,
    random_legal_completer,
    scripted_completer,
)
from llmarena.belief import BeliefTracker
from llmarena.cache import PromptCache
from llmarena.probes import (
    PROBE_KINDS,
    PROBE_SCAFFOLDS,
    Probe,
    _belief_probe,
    _flag_probe,
    _legal_probe,
    balanced,
    generate,
    parse_kinds,
    parse_squares,
    read_jsonl,
    score,
    write_jsonl,
)
from llmarena.view import (
    _grid,
    BELIEF_BLOCK_MARKER,
    FLAG_BLOCK_MARKER,
    LEGAL_BLOCK_MARKER,
    SCAFFOLDS,
    Observation,
    build_observation,
    render,
)

WEIGHTS = BotWeights()


def _permute_hidden_ranks(game: Game, side: Owner, rng: random.Random) -> Game:
    """Relabel every unrevealed piece of ``side``, keeping squares fixed.

    The result is usually an illegal army -- mines outside the rear rows, two
    commanders -- and that is deliberate. Nothing a player may legally see is
    allowed to depend on which rank sits where.
    """
    twin = game.clone()
    squares = [
        square
        for square, piece in twin.board.items()
        if piece.owner == side and not piece.revealed
    ]
    kinds = [twin.board[square].kind for square in squares]
    rng.shuffle(kinds)
    for square, kind in zip(squares, kinds, strict=True):
        twin.board[square] = Piece(side, kind, False)
    return twin


def _play(seed: int, plies: int):
    """Yield ``(game, mover, belief)`` before each of the first ``plies`` moves."""
    game = Game(board=make_opening(seed), turn=Owner(seed % 2))
    players = {
        owner: HeuristicBot(WEIGHTS, seed=seed * 2 + int(owner)) for owner in Owner
    }
    trackers = {owner: BeliefTracker(owner) for owner in Owner}
    for _ in range(plies):
        if game.over:
            return
        mover = game.turn
        yield game, mover, trackers[mover].update(game)
        game.apply(players[mover].choose_move(game))


class ViewTests(unittest.TestCase):
    def test_the_prompt_is_invariant_under_permuting_hidden_ranks(self) -> None:
        rng = random.Random(20260802)
        checked = 0
        for seed in range(4):
            for game, mover, belief in _play(seed, 30):
                baseline = {
                    name: render(build_observation(game, mover, belief), scaffold)
                    for name, scaffold in SCAFFOLDS.items()
                }
                for _ in range(3):
                    twin = _permute_hidden_ranks(game, mover.other, rng)
                    for name, scaffold in SCAFFOLDS.items():
                        self.assertEqual(
                            baseline[name],
                            render(build_observation(twin, mover, belief), scaffold),
                            f"scaffold {name!r} leaked a hidden rank on seed {seed}",
                        )
                checked += 1
        self.assertGreater(checked, 50)

    def test_an_enemy_entry_has_no_rank_field_at_all(self) -> None:
        for game, mover, belief in _play(1, 6):
            observation = build_observation(game, mover, belief)
            for entry in observation.enemy:
                self.assertFalse(hasattr(entry, "kind"))
            payload = observation.to_dict()
            for entry in payload["enemy"]:
                self.assertEqual(set(entry), {"square", "revealed"})

    def test_a_revealed_enemy_piece_is_drawn_as_a_flag(self) -> None:
        # Only a flag is ever revealed, so this is the single enemy rank a
        # prompt is allowed to name.
        game = Game(board=make_opening(5), turn=Owner.HUMAN)
        target = game.flag_candidates(Owner.BOT)[0]
        game.board[target] = Piece(Owner.BOT, PieceKind.FLAG, True)
        observation = build_observation(game, Owner.HUMAN)
        self.assertIn("[旗]", render(observation, SCAFFOLDS["raw"]))
        # Count on the grid, not the whole prompt -- the legend explains the
        # glyph and would otherwise be counted as a second flag.
        self.assertEqual(_grid(observation).count("[旗]"), 1)

    def test_the_observation_survives_a_json_round_trip(self) -> None:
        for game, mover, belief in _play(2, 12):
            observation = build_observation(game, mover, belief)
            restored = Observation.from_dict(
                json.loads(json.dumps(observation.to_dict()))
            )
            self.assertEqual(observation, restored)

    def test_the_scaffolds_are_nested_supersets(self) -> None:
        # A ladder only measures scaffolding if each rung strictly adds text.
        for game, mover, belief in _play(3, 8):
            observation = build_observation(game, mover, belief)
            raw = render(observation, SCAFFOLDS["raw"])
            legal = render(observation, SCAFFOLDS["legal"])
            derived = render(observation, SCAFFOLDS["derived"])
            self.assertNotIn(LEGAL_BLOCK_MARKER, raw)
            self.assertIn(LEGAL_BLOCK_MARKER, legal)
            self.assertIn(FLAG_BLOCK_MARKER, derived)
            self.assertLess(len(raw), len(legal))
            self.assertLessEqual(len(legal), len(derived))

    def test_the_rules_block_is_byte_identical_across_calls(self) -> None:
        # Prompt caching only pays if the shared prefix never varies.
        prefixes = set()
        for game, mover, belief in _play(4, 10):
            text = render(build_observation(game, mover, belief), SCAFFOLDS["legal"])
            prefixes.add(text.split("\n\n")[0])
        self.assertEqual(len(prefixes), 1)


class BeliefTests(unittest.TestCase):
    def test_the_tracker_agrees_with_the_search_bot(self) -> None:
        """Guards the deliberate duplication of the engine's deduction."""
        game = Game(board=make_opening(7), turn=Owner.HUMAN)
        subject = SearchBot(WEIGHTS, seed=1, samples=1, beam_width=3, reply_width=1)
        other = HeuristicBot(WEIGHTS, seed=2)
        tracker = BeliefTracker(Owner.HUMAN)
        compared = 0
        for _ in range(40):
            if game.over:
                break
            if game.turn == Owner.HUMAN:
                move = subject.choose_move(game)
                believed = tracker.update(game)
                self.assertIsNotNone(subject.knowledge)
                self.assertEqual(dict(subject.knowledge.possible), dict(believed))
                compared += 1
                game.apply(move)
            else:
                game.apply(other.choose_move(game))
        self.assertGreater(compared, 10)


class CacheTests(unittest.TestCase):
    def test_a_stored_response_replays(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            cache = PromptCache(Path(folder))
            self.assertIsNone(cache.get("m", "legal", "prompt"))
            cache.put("m", "legal", "prompt", "A10-A9")
            self.assertEqual(cache.get("m", "legal", "prompt"), "A10-A9")
            self.assertEqual(cache.stats()["hits"], 1)

    def test_a_retry_does_not_replay_the_rejected_answer(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            cache = PromptCache(Path(folder))
            cache.put("m", "legal", "prompt", "bad", variant=0)
            self.assertIsNone(cache.get("m", "legal", "prompt", variant=1))

    def test_the_key_separates_model_and_scaffold(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            cache = PromptCache(Path(folder))
            cache.put("a", "legal", "p", "one")
            cache.put("b", "legal", "p", "two")
            cache.put("a", "raw", "p", "three")
            self.assertEqual(cache.get("a", "legal", "p"), "one")
            self.assertEqual(cache.get("b", "legal", "p"), "two")
            self.assertEqual(cache.get("a", "raw", "p"), "three")


class AgentTests(unittest.TestCase):
    def test_a_move_is_parsed_out_of_surrounding_reasoning(self) -> None:
        self.assertEqual(str(parse_response("A10-A9")), "A10-A9")
        self.assertEqual(str(parse_response("先考虑 B7-B8，但我选 C6 -> C7")), "C6-C7")
        self.assertEqual(str(parse_response("答案：\ne2 e3")), "E2-E3")
        self.assertIsNone(parse_response("我不知道"))
        self.assertIsNone(parse_response(""))

    def test_the_agent_always_returns_a_legal_move(self) -> None:
        game = Game(board=make_opening(11), turn=Owner.HUMAN)
        agent = LLMAgent(
            scripted_completer(["完全是胡说八道"]),
            model="stub",
            scaffold=SCAFFOLDS["legal"],
            seed=5,
        )
        for _ in range(8):
            if game.over:
                break
            move = agent.choose_move(game)
            self.assertIn(move, game.legal_moves(game.turn))
            game.apply(move)
        self.assertTrue(all(not r.proposal_legal for r in agent.transcript))
        self.assertTrue(all(r.source == "random" for r in agent.transcript))

    def test_an_illegal_proposal_is_recorded_rather_than_hidden(self) -> None:
        game = Game(board=make_opening(12), turn=Owner.HUMAN)
        # A1-A2 is a real coordinate pair but never legal for the south side.
        agent = LLMAgent(
            scripted_completer(["A1-A2"]),
            model="stub",
            scaffold=SCAFFOLDS["legal"],
            seed=6,
        )
        move = agent.choose_move(game)
        record = agent.transcript[0]
        self.assertEqual(record.proposed, "A1-A2")
        self.assertFalse(record.proposal_legal)
        self.assertEqual(record.source, "random")
        self.assertEqual(record.played, str(move))
        self.assertIn(move, game.legal_moves(Owner.HUMAN))

    def test_forfeit_is_recoverable_from_a_fallback_transcript(self) -> None:
        """One run of paid calls, two penalty regimes."""
        game = Game(board=make_opening(13), turn=Owner.HUMAN)
        legal_first = str(game.legal_moves(Owner.HUMAN)[0])
        agent = LLMAgent(
            scripted_completer([legal_first, "A1-A2"]),
            model="stub",
            scaffold=SCAFFOLDS["legal"],
            seed=7,
        )
        game.apply(agent.choose_move(game))
        game.apply(HeuristicBot(WEIGHTS, seed=1).choose_move(game))
        agent.choose_move(game)
        self.assertIsNone(first_illegal_ply(agent.transcript[:1]))
        self.assertEqual(first_illegal_ply(agent.transcript), 2)

    def test_the_control_completer_never_proposes_an_illegal_move(self) -> None:
        game = Game(board=make_opening(14), turn=Owner.HUMAN)
        agent = LLMAgent(
            random_legal_completer(seed=3),
            model="stub",
            scaffold=SCAFFOLDS["legal"],
            seed=8,
        )
        other = HeuristicBot(WEIGHTS, seed=9)
        for _ in range(20):
            if game.over:
                break
            game.apply(
                agent.choose_move(game)
                if game.turn == Owner.HUMAN
                else other.choose_move(game)
            )
        self.assertTrue(agent.transcript)
        self.assertTrue(all(r.proposal_legal for r in agent.transcript))
        self.assertTrue(all(r.source == "model" for r in agent.transcript))

    def test_two_settings_do_not_share_cache_entries(self) -> None:
        """A re-run at a new effort must not replay the old effort's answers.

        The cache namespace has to carry every knob that changes the reply, not
        just the scaffold -- otherwise a sweep over `effort` would return the
        first setting's moves for all of them and print a spotless null result.
        """
        with tempfile.TemporaryDirectory() as folder:
            cache = PromptCache(Path(folder))
            game = Game(board=make_opening(16), turn=Owner.HUMAN)
            legal = game.legal_moves(Owner.HUMAN)
            first, second = str(legal[0]), str(legal[1])

            def agent_for(namespace: str, reply: str) -> LLMAgent:
                return LLMAgent(
                    scripted_completer([reply]),
                    model="stub",
                    scaffold=SCAFFOLDS["legal"],
                    cache=cache,
                    cache_namespace=namespace,
                    seed=1,
                )

            low = agent_for("legal/low/adaptive", first)
            self.assertEqual(str(low.choose_move(game)), first)
            high = agent_for("legal/high/adaptive", second)
            self.assertEqual(str(high.choose_move(game)), second)
            # Same namespace replays rather than re-asking.
            again = agent_for("legal/low/adaptive", second)
            self.assertEqual(str(again.choose_move(game)), first)

    def test_the_cache_makes_a_replay_free(self) -> None:
        with tempfile.TemporaryDirectory() as folder:
            cache = PromptCache(Path(folder))
            calls = {"n": 0}

            def counting(prompt: str) -> str:
                calls["n"] += 1
                return random_legal_completer(seed=1)(prompt)

            def run() -> list[str]:
                game = Game(board=make_opening(15), turn=Owner.HUMAN)
                agent = LLMAgent(
                    counting,
                    model="stub",
                    scaffold=SCAFFOLDS["legal"],
                    cache=cache,
                    seed=2,
                )
                other = HeuristicBot(WEIGHTS, seed=4)
                played = []
                for _ in range(10):
                    if game.over:
                        break
                    if game.turn == Owner.HUMAN:
                        move = agent.choose_move(game)
                        played.append(str(move))
                    else:
                        move = other.choose_move(game)
                    game.apply(move)
                return played

            first = run()
            spent = calls["n"]
            self.assertGreater(spent, 0)
            self.assertEqual(run(), first)
            self.assertEqual(calls["n"], spent)  # second run paid nothing


class ProbeTests(unittest.TestCase):
    def test_a_legal_move_probe_label_is_the_engines_destination_set(self) -> None:
        rng = random.Random(1)
        found = 0
        for game, mover, belief in _play(21, 30):
            observation = build_observation(game, mover, belief)
            probe = _legal_probe(observation, rng, "x")
            if probe is None:
                continue
            square = next(
                p.square
                for p in observation.own
                if format_position(p.square) == probe.meta["square"]
            )
            truth = {
                format_position(target)
                for target in game._destinations(square, set(game.board))
            }
            self.assertEqual(set(probe.label), truth)
            found += 1
        self.assertGreater(found, 5)

    def test_a_flag_probe_label_is_the_engines_candidate_set(self) -> None:
        found = 0
        for game, mover, belief in _play(22, 30):
            observation = build_observation(game, mover, belief)
            probe = _flag_probe(observation, "x")
            self.assertIsNotNone(probe)
            assert probe is not None
            truth = {
                format_position(s) for s in game.flag_candidates(mover.other)
            }
            self.assertEqual(set(probe.label), truth)
            found += 1
        self.assertGreater(found, 5)

    def test_a_belief_probe_label_is_the_engines_deduced_set(self) -> None:
        rng = random.Random(2)
        found = 0
        for seed in range(30):
            for game, mover, belief in _play(seed, 60):
                observation = build_observation(game, mover, belief)
                probe = _belief_probe(observation, rng, "x")
                if probe is None:
                    continue
                square = next(
                    s
                    for s in belief
                    if format_position(s) == probe.meta["square"]
                )
                self.assertEqual(
                    set(probe.label), {k.name for k in belief[square]}
                )
                found += 1
            if found > 5:
                break
        self.assertGreater(found, 5)

    def test_the_true_rank_is_recorded_for_scoring_but_never_rendered(self) -> None:
        """`meta` must not reach the model.

        The belief label is the engine's per-square deduction, which does no
        rank counting -- so it is a loose upper bound and a better answer can
        be a strict subset of it. Scoring needs the piece's actual rank; the
        prompt must not have it.
        """
        probes = [
            p
            for p in generate(range(61, 75), WEIGHTS, max_plies=60)
            if p.kind == "belief"
        ]
        self.assertTrue(probes)
        for probe in probes:
            self.assertIn("true_kind", probe.meta)
            self.assertIsNotNone(probe.meta["true_kind"])
            # The true rank is always inside the engine's set: the engine is
            # loose, never wrong.
            self.assertIn(probe.meta["true_kind"], probe.label)
            # prompt() is a pure function of observation + question.
            twin = Probe(
                probe_id=probe.probe_id,
                kind=probe.kind,
                observation=probe.observation,
                question=probe.question,
                label=probe.label,
                meta={"square": probe.meta["square"], "true_kind": "COMMANDER"},
            )
            self.assertEqual(probe.prompt(), twin.prompt())

    def test_out_deducing_the_engine_scores_as_correct_not_as_an_error(self) -> None:
        probe = Probe(
            probe_id="x",
            kind="belief",
            observation=build_observation(
                Game(board=make_opening(42), turn=Owner.HUMAN), Owner.HUMAN
            ),
            question="q",
            label=("BRIGADIER", "GENERAL", "MAJOR_GENERAL"),
            meta={"square": "D12", "true_kind": "MAJOR_GENERAL"},
        )
        # Drops GENERAL, keeps the truth -- tighter than the engine, still right.
        tighter = score(probe, "师 旅")
        self.assertFalse(tighter["exact"])  # equality scoring calls this wrong
        self.assertTrue(tighter["correct"])  # ground truth calls it right
        self.assertEqual(tighter["tighter_by"], 1)

        # Dropping the true rank is a real error.
        wrong = score(probe, "旅")
        self.assertFalse(wrong["sound"])
        self.assertFalse(wrong["correct"])

        # Adding a rank the engine proved impossible is also a real error.
        wide = score(probe, "师 旅 军 工")
        self.assertTrue(wide["sound"])
        self.assertTrue(wide["over_wide"])
        self.assertFalse(wide["correct"])

    def test_a_probe_prompt_never_contains_its_own_answer(self) -> None:
        forbidden = {
            "legal_moves": LEGAL_BLOCK_MARKER,
            "flag_candidates": FLAG_BLOCK_MARKER,
            "belief": BELIEF_BLOCK_MARKER,
            "tightest": BELIEF_BLOCK_MARKER,
        }
        probes = list(generate(range(24, 30), WEIGHTS, max_plies=50))
        self.assertTrue(probes)
        seen = set()
        for probe in probes:
            seen.add(probe.kind)
            self.assertNotIn(forbidden[probe.kind], probe.prompt())
        self.assertEqual(seen, set(forbidden))

    def test_probe_scaffolds_switch_off_exactly_the_answering_block(self) -> None:
        self.assertFalse(PROBE_SCAFFOLDS["legal_moves"].legal_moves)
        self.assertFalse(PROBE_SCAFFOLDS["flag_candidates"].flag_candidates)
        self.assertFalse(PROBE_SCAFFOLDS["belief"].belief)

    def test_the_battery_survives_a_json_round_trip(self) -> None:
        probes = list(generate(range(31, 34), WEIGHTS, max_plies=40))
        self.assertTrue(probes)
        with tempfile.TemporaryDirectory() as folder:
            path = write_jsonl(probes, Path(folder) / "probes.jsonl")
            restored = read_jsonl(path)
        self.assertEqual(len(restored), len(probes))
        for before, after in zip(probes, restored, strict=True):
            self.assertEqual(before.to_dict(), after.to_dict())
            self.assertEqual(before.prompt(), after.prompt())

    def test_answers_are_parsed_from_the_final_line_only(self) -> None:
        self.assertEqual(
            parse_squares("先想想 A1 B2\n最终答案\nC3 D4"), {"C3", "D4"}
        )
        self.assertEqual(parse_kinds("分析中\n司 军 工"), {"COMMANDER", "GENERAL", "ENGINEER"})

    def test_scoring_separates_exact_from_partial(self) -> None:
        probe = Probe(
            probe_id="x",
            kind="flag_candidates",
            observation=build_observation(
                Game(board=make_opening(41), turn=Owner.HUMAN), Owner.HUMAN
            ),
            question="q",
            label=("B1", "D1"),
            meta={},
        )
        self.assertTrue(score(probe, "B1 D1")["exact"])
        partial = score(probe, "B1")
        self.assertFalse(partial["exact"])
        self.assertAlmostEqual(partial["jaccard"], 0.5)
        self.assertEqual(score(probe, "无")["jaccard"], 0.0)

    def test_the_battery_is_balanced_across_kinds(self) -> None:
        probes = list(generate(range(51, 71), WEIGHTS, max_plies=80))
        chosen = balanced(probes, per_kind=5)
        counts = {kind: 0 for kind in PROBE_KINDS}
        for probe in chosen:
            counts[probe.kind] += 1
        for kind, count in counts.items():
            self.assertGreater(count, 0, f"no {kind} probes were kept")
            self.assertLessEqual(count, 5)


class ArenaIntegrationTests(unittest.TestCase):
    SPEC = AgentSpec(
        "llm-stub",
        "external",
        builder="llmarena.agent:build_agent",
        options=(
            ("completer", "random-legal"),
            ("scaffold", "legal"),
            ("model", "stub"),
        ),
    )

    def test_a_job_carrying_an_external_subject_pickles(self) -> None:
        job = Job(
            WEIGHTS,
            AgentSpec("random", "random"),
            seed=3,
            subject_side=int(Owner.HUMAN),
            subject_spec=self.SPEC,
        )
        self.assertEqual(pickle.loads(pickle.dumps(job)), job)

    def test_the_arena_plays_a_subject_that_is_not_a_search_bot(self) -> None:
        result = play_match(
            Job(
                WEIGHTS,
                AgentSpec("random", "random"),
                seed=3,
                subject_side=int(Owner.HUMAN),
                subject_spec=self.SPEC,
            )
        )
        self.assertIn(result.result, (0.0, 0.5, 1.0))
        self.assertGreater(result.plies, 0)

    def test_omitting_the_subject_spec_still_builds_a_search_bot(self) -> None:
        # The default path is what every historical measurement used; it must
        # not change shape just because the seat became configurable.
        result = play_match(
            Job(
                WEIGHTS,
                AgentSpec("random", "random"),
                seed=3,
                subject_side=int(Owner.HUMAN),
            )
        )
        self.assertIn(result.result, (0.0, 0.5, 1.0))


if __name__ == "__main__":
    unittest.main()
