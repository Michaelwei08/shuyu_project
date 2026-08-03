"""The real completer: one Junqi position in, one move out.

``anthropic`` is imported lazily inside :func:`_client` so the rest of the
package -- the view, the probes, the cache, every test -- keeps working with no
SDK installed and no credentials. Only this function needs either.

Three API facts drive the shape of this file, and each one would be a silent
bug if assumed rather than looked up:

* **Sampling parameters are gone.** ``temperature`` / ``top_p`` / ``top_k``
  return a 400 on Opus 5 and Sonnet 5. There is therefore *no* way to ask for a
  deterministic sample, which is why :mod:`llmarena.cache` is not a
  nice-to-have: a warm cache is the only reproducibility this harness gets.
* **Thinking is on by default on Opus 5**, and ``max_tokens`` caps thinking
  *plus* answer. A tight ``max_tokens`` truncates the move mid-sentence. It is
  also gated: ``thinking={"type": "disabled"}`` is only legal at effort
  ``high`` or below, so that combination is rejected here rather than at the
  API.
* **The cached prefix has a per-model minimum** -- 512 tokens on Opus 5, 1024
  on Sonnet 5, 4096 on Haiku 4.5 -- and a prefix under it silently does not
  cache. ``RULES`` is ~600-800 tokens, so it caches on Opus 5 and very likely
  does not on the cheaper models. Nothing errors either way, so every response's
  ``cache_read_input_tokens`` is recorded and the runners print it.

Prices used by :mod:`llmarena.cost` were current 2026-06-24; re-check before
quoting a bill.
"""

from __future__ import annotations

import threading
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any

from .view import RULES

#: Effort levels that forbid disabling thinking (the API returns 400).
_EFFORT_FORBIDS_DISABLED = frozenset({"xhigh", "max"})
_EFFORTS = ("low", "medium", "high", "xhigh", "max")

#: Beta flag for the `fallbacks: "default"` scalar form. The array form uses a
#: different, earlier header; pairing either header with the other shape 400s.
_FALLBACK_BETA = "server-side-fallback-2026-07-01"


@dataclass
class Usage:
    """Token totals across a run, aggregated across threads."""

    calls: int = 0
    input_tokens: int = 0
    output_tokens: int = 0
    cache_read_tokens: int = 0
    cache_write_tokens: int = 0
    refusals: int = 0
    _lock: threading.Lock = field(default_factory=threading.Lock, repr=False)

    def record(self, usage: Any, refused: bool = False) -> None:
        with self._lock:
            self.calls += 1
            self.input_tokens += getattr(usage, "input_tokens", 0) or 0
            self.output_tokens += getattr(usage, "output_tokens", 0) or 0
            self.cache_read_tokens += (
                getattr(usage, "cache_read_input_tokens", 0) or 0
            )
            self.cache_write_tokens += (
                getattr(usage, "cache_creation_input_tokens", 0) or 0
            )
            if refused:
                self.refusals += 1

    def format(self) -> str:
        cached = self.cache_read_tokens
        total_in = self.input_tokens + cached + self.cache_write_tokens
        share = f"{cached / total_in:.0%}" if total_in else "n/a"
        return (
            f"{self.calls} calls | in {total_in:,} "
            f"(uncached {self.input_tokens:,}, cache-read {cached:,} = {share}, "
            f"cache-write {self.cache_write_tokens:,}) | "
            f"out {self.output_tokens:,} | refusals {self.refusals}"
        )


def split_cacheable(prompt: str) -> tuple[str | None, str]:
    """Split a rendered prompt into its cacheable prefix and the rest.

    ``render`` puts ``RULES`` first and joins blocks with a blank line, so the
    prefix is recognised by an exact prefix match rather than a sentinel. A
    prompt rendered with ``Scaffold(rules=False)`` simply has no cacheable part.
    """
    if prompt.startswith(RULES):
        return RULES, prompt[len(RULES) :].lstrip("\n")
    return None, prompt


_CLIENTS: dict[tuple[str, int], Any] = {}
_CLIENTS_LOCK = threading.Lock()


def _client(max_retries: int) -> Any:
    """One SDK client per thread-safe config, built on first use."""
    key = ("default", max_retries)
    with _CLIENTS_LOCK:
        if key not in _CLIENTS:
            try:
                import anthropic
            except ModuleNotFoundError as error:  # pragma: no cover - env dependent
                raise ModuleNotFoundError(
                    "the anthropic SDK is required to call a real model: "
                    "pip install anthropic"
                ) from error
            # A rate-limit storm must not surface as a failed match -- the arena
            # aborts a whole run above a 2% failure rate. The SDK's own backoff
            # is the right place to absorb it.
            _CLIENTS[key] = anthropic.Anthropic(max_retries=max_retries)
        return _CLIENTS[key]


def anthropic_completer(
    model: str = "claude-opus-5",
    *,
    effort: str = "low",
    thinking: str = "adaptive",
    max_tokens: int = 8000,
    max_retries: int = 8,
    fallbacks: bool = True,
    usage: Usage | None = None,
) -> Callable[[str], str]:
    """Build a ``(prompt) -> text`` completer backed by the Messages API.

    ``effort`` is the main cost lever and a legitimate experimental axis in its
    own right: ``low`` on this model is strong, and it is a far better knob than
    dropping to a weaker model, which would also lose prompt caching.

    ``thinking`` is ``"adaptive"`` or ``"disabled"``. Disabling it is only legal
    at effort ``high`` or below, and on Opus 5 it can leak ``<thinking>`` tags
    into the visible text -- harmless here, because :func:`parse_response` reads
    the *last* coordinate pair, but the instruction below discourages it anyway.
    """
    if effort not in _EFFORTS:
        raise ValueError(f"effort must be one of {_EFFORTS}, got {effort!r}")
    if thinking not in ("adaptive", "disabled"):
        raise ValueError("thinking must be 'adaptive' or 'disabled'")
    if thinking == "disabled" and effort in _EFFORT_FORBIDS_DISABLED:
        raise ValueError(
            f"thinking='disabled' is rejected at effort={effort!r}; "
            "use effort 'high' or lower, or leave thinking adaptive"
        )

    def complete(prompt: str) -> str:
        cacheable, rest = split_cacheable(prompt)
        system: list[dict[str, Any]] = []
        if cacheable is not None:
            system.append(
                {
                    "type": "text",
                    "text": cacheable,
                    # The one stable prefix in the whole request. Whether it
                    # actually caches depends on the model's minimum -- read
                    # `cache_read_input_tokens` rather than assuming.
                    "cache_control": {"type": "ephemeral"},
                }
            )

        request: dict[str, Any] = {
            "model": model,
            "max_tokens": max_tokens,
            "messages": [{"role": "user", "content": rest}],
            "thinking": {"type": thinking},
            "output_config": {"effort": effort},
        }
        if system:
            request["system"] = system
        # No temperature / top_p / top_k -- removed on this model family, 400.

        client = _client(max_retries)
        if fallbacks:
            response = client.beta.messages.create(
                betas=[_FALLBACK_BETA], fallbacks="default", **request
            )
        else:
            response = client.messages.create(**request)

        refused = response.stop_reason == "refusal"
        if usage is not None:
            usage.record(response.usage, refused=refused)
        if refused:
            # A declined request returns HTTP 200 with empty or partial content.
            # Returning "" makes the agent fall back to a legal move and records
            # the ply as an unparseable proposal, which is the honest scoring.
            return ""
        return "".join(
            block.text for block in response.content if block.type == "text"
        )

    return complete


#: Shared by every completer built through the ``AgentSpec`` registry. A match
#: run builds one agent per game, so without a shared accumulator each would
#: report its own handful of calls and the run total would be lost.
USAGE = Usage()


def register() -> None:
    """Add ``anthropic`` to the completer registry used by ``AgentSpec``."""
    from .agent import COMPLETERS

    def build(options: dict[str, str]):
        return anthropic_completer(
            model=options.get("model", "claude-opus-5"),
            effort=options.get("effort", "low"),
            thinking=options.get("thinking", "adaptive"),
            max_tokens=int(options.get("max_tokens", "8000")),
            fallbacks=options.get("fallbacks", "1") not in ("0", "false", "no"),
            usage=USAGE,
        )

    COMPLETERS["anthropic"] = build
