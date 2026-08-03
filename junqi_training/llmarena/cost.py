"""Price a run before paying for it.

    python -m llmarena.cost --probes data/probes.jsonl --games 60

Token counts here are a **character-based estimate**, not a measurement: the
real number comes from ``client.messages.count_tokens``, which needs
credentials this module deliberately does not require. Treat the output as an
order of magnitude for deciding whether to proceed, and re-measure before
quoting anything.

The per-model cache minimum is the load-bearing column. A prefix under it does
not cache, silently -- no error, ``cache_creation_input_tokens: 0`` -- so the
same code is ~30% cheaper per call on one model and full price on another.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from .view import RULES

#: Per-million-token prices, current 2026-06-24. Re-check before quoting.
#: ``cache_min`` is the smallest prefix that will cache on that model.
@dataclass(frozen=True)
class Model:
    name: str
    input_per_mtok: float
    output_per_mtok: float
    cache_min_tokens: int


MODELS: dict[str, Model] = {
    "claude-opus-5": Model("claude-opus-5", 5.00, 25.00, 512),
    "claude-opus-4-8": Model("claude-opus-4-8", 5.00, 25.00, 1024),
    # Introductory pricing runs through 2026-08-31; standard is 3.00 / 15.00.
    "claude-sonnet-5": Model("claude-sonnet-5", 2.00, 10.00, 1024),
    "claude-haiku-4-5": Model("claude-haiku-4-5", 1.00, 5.00, 4096),
    "claude-fable-5": Model("claude-fable-5", 10.00, 50.00, 512),
}

CACHE_READ_MULTIPLIER = 0.1
CACHE_WRITE_MULTIPLIER = 1.25
BATCH_MULTIPLIER = 0.5

#: Codepoint ranges counted as roughly one token per character.
_WIDE_RANGES = (
    (0x3000, 0x303F),  # CJK punctuation
    (0x4E00, 0x9FFF),  # CJK unified ideographs
    (0xFF00, 0xFFEF),  # fullwidth forms
)


def estimate_tokens(text: str) -> int:
    """Rough token count: ~1 per CJK character, ~1 per 3.5 ASCII characters.

    Deliberately crude. Its only job is to answer "is this run ten dollars or a
    thousand" without an API key.
    """
    wide = sum(
        1
        for character in text
        if any(low <= ord(character) <= high for low, high in _WIDE_RANGES)
    )
    return round(wide + (len(text) - wide) / 3.5)


@dataclass
class Estimate:
    model: str
    calls: int
    input_tokens: int
    cached_tokens: int
    output_tokens: int
    dollars: float
    caches: bool

    def format(self) -> str:
        note = "" if self.caches else "  (prefix below cache minimum -- no caching)"
        return (
            f"  {self.model:<18} {self.calls:>7,} calls  "
            f"in {self.input_tokens + self.cached_tokens:>10,}  "
            f"out {self.output_tokens:>8,}  ${self.dollars:>8.2f}{note}"
        )


def estimate(
    model: str,
    calls: int,
    prompt_tokens: int,
    output_tokens: int,
    cacheable_tokens: int = 0,
    batch: bool = False,
) -> Estimate:
    """Cost of ``calls`` identical-shaped requests.

    Assumes the cacheable prefix is written once and read on every later call,
    which is the best case; a cache that expires between calls costs more.
    """
    spec = MODELS[model]
    caches = cacheable_tokens >= spec.cache_min_tokens
    cached = cacheable_tokens if caches else 0
    uncached_per_call = prompt_tokens - cached

    input_cost = uncached_per_call * calls
    if caches:
        input_cost += cached * CACHE_WRITE_MULTIPLIER  # one write
        input_cost += cached * CACHE_READ_MULTIPLIER * max(0, calls - 1)

    dollars = (
        input_cost / 1_000_000 * spec.input_per_mtok
        + output_tokens * calls / 1_000_000 * spec.output_per_mtok
    )
    if batch:
        dollars *= BATCH_MULTIPLIER
    return Estimate(
        model=model,
        calls=calls,
        input_tokens=uncached_per_call * calls,
        cached_tokens=cached * calls,
        output_tokens=output_tokens * calls,
        dollars=dollars,
        caches=caches,
    )


#: Mean plies per game measured over 12 heuristic-vs-heuristic games; the LLM
#: moves on about half of them. The distribution has a heavy tail -- 1 game in
#: 12 hit the 300-ply cap -- so budget with headroom.
MEAN_PLIES = 85
LLM_SHARE = 0.5


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probes", type=Path, default=Path("data/probes.jsonl"))
    parser.add_argument("--scaffolds", type=int, default=3)
    parser.add_argument("--games", type=int, default=60)
    parser.add_argument("--match-scaffolds", type=int, default=2)
    parser.add_argument(
        "--output-tokens",
        type=int,
        default=600,
        help="mean output per call INCLUDING thinking, which is billed as output",
    )
    parser.add_argument("--batch-probes", action="store_true")
    args = parser.parse_args()

    from .probes import read_jsonl

    rules_tokens = estimate_tokens(RULES)
    if args.probes.exists():
        probes = read_jsonl(args.probes)
        prompt_tokens = round(
            sum(estimate_tokens(p.prompt()) for p in probes) / len(probes)
        )
        probe_calls = len(probes) * args.scaffolds
    else:
        print(f"note: {args.probes} not found; using a 1500-token placeholder\n")
        prompt_tokens, probe_calls = 1500, 300 * args.scaffolds

    match_calls = round(args.games * MEAN_PLIES * LLM_SHARE) * args.match_scaffolds

    print(f"cacheable prefix (RULES): ~{rules_tokens} tokens")
    print(f"mean prompt:              ~{prompt_tokens} tokens")
    print(f"assumed output per call:   {args.output_tokens} tokens (incl. thinking)")
    print("\nestimates are character-based, not measured -- treat as +/- 25%\n")

    for label, calls, batch in (
        (f"probe battery x{args.scaffolds} scaffolds", probe_calls, args.batch_probes),
        (f"{args.games} games x{args.match_scaffolds} scaffolds", match_calls, False),
    ):
        print(f"{label}{'  [batch, 50% off]' if batch else ''}:")
        for name in MODELS:
            print(
                estimate(
                    name,
                    calls,
                    prompt_tokens,
                    args.output_tokens,
                    cacheable_tokens=rules_tokens,
                    batch=batch,
                ).format()
            )
        print()


if __name__ == "__main__":
    main()
