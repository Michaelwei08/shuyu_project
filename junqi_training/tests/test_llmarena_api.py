"""Tests for the paid path, exercised against a fake client.

Nothing here touches the network or needs the `anthropic` SDK installed. The
point is to pin the three request-shape facts that would each be a silent 400
or a silent overcharge in production: no sampling parameters, the rules block
sent as a cached system block, and disabled-thinking rejected above `high`
effort before the API sees it.
"""

import types
import unittest

from llmarena import anthropic_completer as api
from llmarena.cost import MODELS, estimate, estimate_tokens
from llmarena.run_matches import rollouts, search_equivalent
from llmarena.view import RULES, SCAFFOLDS, build_observation, render


def _response(text="A10-A9", stop_reason="end_turn", **usage):
    totals = {
        "input_tokens": 100,
        "output_tokens": 20,
        "cache_read_input_tokens": 0,
        "cache_creation_input_tokens": 0,
    }
    totals.update(usage)
    return types.SimpleNamespace(
        stop_reason=stop_reason,
        content=[types.SimpleNamespace(type="text", text=text)],
        usage=types.SimpleNamespace(**totals),
    )


class FakeClient:
    """Records request bodies and replays a canned response."""

    def __init__(self, response=None):
        self.requests: list[dict] = []
        self.beta_requests: list[dict] = []
        self.response = response or _response()

        def create(**kwargs):
            self.requests.append(kwargs)
            return self.response

        def beta_create(**kwargs):
            self.beta_requests.append(kwargs)
            return self.response

        self.messages = types.SimpleNamespace(create=create)
        self.beta = types.SimpleNamespace(
            messages=types.SimpleNamespace(create=beta_create)
        )


class CompleterTests(unittest.TestCase):
    def setUp(self) -> None:
        self._real_client = api._client
        self.client = FakeClient()
        api._client = lambda _retries: self.client

    def tearDown(self) -> None:
        api._client = self._real_client

    def _prompt(self) -> str:
        from junqi.arena import make_opening
        from junqi.game import Game
        from junqi.types import Owner

        game = Game(board=make_opening(1), turn=Owner.HUMAN)
        return render(build_observation(game, Owner.HUMAN), SCAFFOLDS["legal"])

    def test_the_cacheable_prefix_is_split_off(self) -> None:
        prefix, rest = api.split_cacheable(self._prompt())
        self.assertEqual(prefix, RULES)
        self.assertNotIn(RULES, rest)
        self.assertTrue(rest.startswith("你执"))

    def test_a_prompt_without_the_rules_has_no_cacheable_prefix(self) -> None:
        from junqi.arena import make_opening
        from junqi.game import Game
        from junqi.types import Owner
        from llmarena.view import Scaffold

        game = Game(board=make_opening(1), turn=Owner.HUMAN)
        bare = render(
            build_observation(game, Owner.HUMAN), Scaffold("bare", rules=False)
        )
        prefix, rest = api.split_cacheable(bare)
        self.assertIsNone(prefix)
        self.assertEqual(rest, bare)

    def test_the_request_carries_no_sampling_parameters(self) -> None:
        # temperature / top_p / top_k are removed on this model family and
        # return a 400. There is no deterministic-sampling option to add back.
        api.anthropic_completer(usage=None)(self._prompt())
        body = self.client.beta_requests[0]
        for banned in ("temperature", "top_p", "top_k"):
            self.assertNotIn(banned, body)

    def test_the_rules_block_is_sent_as_a_cached_system_block(self) -> None:
        api.anthropic_completer()(self._prompt())
        body = self.client.beta_requests[0]
        self.assertEqual(len(body["system"]), 1)
        block = body["system"][0]
        self.assertEqual(block["text"], RULES)
        self.assertEqual(block["cache_control"], {"type": "ephemeral"})
        # The volatile position must be in the user turn, after the breakpoint.
        self.assertNotIn(RULES, body["messages"][0]["content"])

    def test_effort_and_thinking_reach_the_request(self) -> None:
        api.anthropic_completer(effort="medium", thinking="disabled")(self._prompt())
        body = self.client.beta_requests[0]
        self.assertEqual(body["output_config"], {"effort": "medium"})
        self.assertEqual(body["thinking"], {"type": "disabled"})

    def test_disabling_thinking_above_high_effort_is_rejected_locally(self) -> None:
        # The API returns 400 for this pair; catching it here turns a failed
        # run into a failed argument parse.
        for effort in ("xhigh", "max"):
            with self.assertRaises(ValueError):
                api.anthropic_completer(effort=effort, thinking="disabled")
        api.anthropic_completer(effort="high", thinking="disabled")  # legal

    def test_an_unknown_effort_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            api.anthropic_completer(effort="extreme")

    def test_fallbacks_choose_the_beta_endpoint(self) -> None:
        api.anthropic_completer(fallbacks=True)(self._prompt())
        self.assertEqual(len(self.client.beta_requests), 1)
        self.assertEqual(len(self.client.requests), 0)

        api.anthropic_completer(fallbacks=False)(self._prompt())
        self.assertEqual(len(self.client.requests), 1)

    def test_a_refusal_returns_empty_text_rather_than_raising(self) -> None:
        # HTTP 200 with stop_reason "refusal" and empty content; reading
        # content[0] unconditionally would raise. The agent then falls back to
        # a legal move and the ply is recorded as an unparseable proposal.
        self.client.response = types.SimpleNamespace(
            stop_reason="refusal",
            content=[],
            usage=types.SimpleNamespace(
                input_tokens=10,
                output_tokens=0,
                cache_read_input_tokens=0,
                cache_creation_input_tokens=0,
            ),
        )
        usage = api.Usage()
        self.assertEqual(api.anthropic_completer(usage=usage)(self._prompt()), "")
        self.assertEqual(usage.refusals, 1)

    def test_usage_accumulates_across_calls(self) -> None:
        self.client.response = _response(
            input_tokens=200, output_tokens=50, cache_read_input_tokens=600
        )
        usage = api.Usage()
        complete = api.anthropic_completer(usage=usage)
        for _ in range(3):
            complete(self._prompt())
        self.assertEqual(usage.calls, 3)
        self.assertEqual(usage.input_tokens, 600)
        self.assertEqual(usage.output_tokens, 150)
        self.assertEqual(usage.cache_read_tokens, 1800)
        self.assertIn("3 calls", usage.format())

    def test_the_agent_drives_the_completer_end_to_end(self) -> None:
        from junqi.arena import make_opening
        from junqi.game import Game
        from junqi.types import Owner
        from llmarena.agent import LLMAgent

        game = Game(board=make_opening(2), turn=Owner.HUMAN)
        legal = game.legal_moves(Owner.HUMAN)[0]
        self.client.response = _response(text=f"我走 {legal}")
        agent = LLMAgent(
            api.anthropic_completer(),
            model="fake",
            scaffold=SCAFFOLDS["legal"],
            seed=1,
        )
        self.assertEqual(agent.choose_move(game), legal)
        self.assertTrue(agent.transcript[0].proposal_legal)


class RegistryTests(unittest.TestCase):
    def test_the_anthropic_backend_registers_on_demand(self) -> None:
        from llmarena.agent import COMPLETERS, make_completer

        COMPLETERS.pop("anthropic", None)
        make_completer("anthropic", {"model": "claude-opus-5"})
        self.assertIn("anthropic", COMPLETERS)

    def test_an_unknown_completer_still_raises(self) -> None:
        from llmarena.agent import make_completer

        with self.assertRaises(ValueError):
            make_completer("gpt-9", {})


class CostTests(unittest.TestCase):
    def test_chinese_costs_far_more_tokens_per_character_than_ascii(self) -> None:
        self.assertGreater(estimate_tokens("军棋" * 100), estimate_tokens("ab" * 100))
        self.assertGreater(estimate_tokens(RULES), 400)

    def test_a_prefix_below_the_model_minimum_does_not_cache(self) -> None:
        # 700 tokens clears Opus 5's 512 minimum and misses Haiku's 4096 --
        # the same code is cheaper on one model and full price on the other.
        opus = estimate("claude-opus-5", 1000, 1500, 500, cacheable_tokens=700)
        haiku = estimate("claude-haiku-4-5", 1000, 1500, 500, cacheable_tokens=700)
        self.assertTrue(opus.caches)
        self.assertFalse(haiku.caches)
        self.assertEqual(haiku.cached_tokens, 0)

    def test_caching_lowers_the_bill(self) -> None:
        without = estimate("claude-opus-5", 500, 1500, 400, cacheable_tokens=0)
        with_cache = estimate("claude-opus-5", 500, 1500, 400, cacheable_tokens=700)
        self.assertLess(with_cache.dollars, without.dollars)

    def test_batch_halves_the_bill(self) -> None:
        plain = estimate("claude-opus-5", 100, 1500, 400)
        batched = estimate("claude-opus-5", 100, 1500, 400, batch=True)
        self.assertAlmostEqual(batched.dollars, plain.dollars / 2)

    def test_every_model_has_a_cache_minimum(self) -> None:
        for spec in MODELS.values():
            self.assertGreater(spec.cache_min_tokens, 0)


class SearchEquivalentTests(unittest.TestCase):
    def test_rollouts_are_monotone_across_the_ladder(self) -> None:
        self.assertLess(rollouts(1, 4, 1), rollouts(3, 8, 3))
        self.assertLess(rollouts(3, 8, 3), rollouts(6, 10, 4))

    def test_a_crossing_is_interpolated_between_rungs(self) -> None:
        answer = search_equivalent([(4, 0.8), (72, 0.3), (240, 0.1)])
        self.assertTrue(answer.startswith("~"))
        value = float(answer.lstrip("~").split()[0])
        self.assertGreater(value, 4)
        self.assertLess(value, 72)

    def test_beating_or_losing_to_everything_is_reported_as_a_bound(self) -> None:
        self.assertIn("> 240", search_equivalent([(4, 0.9), (72, 0.8), (240, 0.7)]))
        self.assertIn("< 4", search_equivalent([(4, 0.2), (72, 0.1), (240, 0.0)]))

    def test_a_non_monotone_ladder_is_reported_rather_than_fitted(self) -> None:
        answer = search_equivalent([(4, 0.3), (72, 0.6), (240, 0.55)])
        self.assertIn("non-monotone", answer)


if __name__ == "__main__":
    unittest.main()
