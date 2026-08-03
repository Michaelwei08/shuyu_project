"""Score a model on the diagnostic battery.

    python -m llmarena.run_probes --model claude-opus-5 --effort low --limit 60

One call per probe, graded against a label the engine computed exactly. This is
the cheap instrument: a full match spends ~42 calls to return one win/loss bit,
while the battery returns a graded answer per call and splits the result into
rule execution / public deduction / private inference.

**There is deliberately no move-regret instrument here.** Scoring an LLM's move
against a strong `SearchBot` referee was the obvious third option, and it is the
one thing this repo has already learned not to trust: `CLAUDE.md` records that
`search-mid` -- the pool's ceiling -- is the bot's own policy class, which is
why the pool could not grade anything past parity-with-self until an oracle was
added. A referee drawn from that same policy class would measure agreement with
the bot, not quality. Scaffolding is therefore ablated on real matches
(`run_matches`), and the battery measures capabilities that have ground truth.

Batch API note: these calls are offline and order-independent, so the Batches
API would halve their cost. It is not implemented, because it cannot be
exercised from this machine (no credentials) and untested billing code is worse
than none. It is the obvious next saving if the battery gets re-run across
several models.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from statistics import fmean

from .anthropic_completer import Usage, anthropic_completer
from .cache import NullCache, PromptCache
from .probes import PROBE_KINDS, read_jsonl, score


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probes", type=Path, default=Path("data/probes.jsonl"))
    parser.add_argument(
        "--backend",
        default="anthropic",
        choices=("anthropic", "claude-cli"),
        help="claude-cli runs on Claude Code's own credentials (no API spend) "
        "but puts ~5.1k tokens of harness between the experiment and the model",
    )
    parser.add_argument(
        "--workdir",
        type=Path,
        default=Path("data/cli-scratch"),
        help="empty directory for the claude-cli backend; it refuses to run "
        "anywhere a CLAUDE.md would be auto-loaded into the prompt",
    )
    parser.add_argument("--model", default="claude-opus-5")
    parser.add_argument("--effort", default="low")
    parser.add_argument("--thinking", default="adaptive", choices=("adaptive", "disabled"))
    parser.add_argument("--max-tokens", type=int, default=8000)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument(
        "--limit", type=int, help="score only the first N probes of each kind"
    )
    parser.add_argument("--cache", type=Path, default=Path("data/cache"))
    parser.add_argument("--no-cache", action="store_true")
    parser.add_argument("--out", type=Path, default=Path("data/probe-results.jsonl"))
    args = parser.parse_args()

    probes = read_jsonl(args.probes)
    if args.limit:
        kept: dict[str, int] = defaultdict(int)
        selected = []
        for probe in probes:
            if kept[probe.kind] < args.limit:
                kept[probe.kind] += 1
                selected.append(probe)
        probes = selected
    if not probes:
        raise SystemExit(f"no probes in {args.probes}")

    usage = Usage()
    if args.backend == "claude-cli":
        from .cli_completer import claude_cli_completer

        complete = claude_cli_completer(
            args.workdir, model=args.model, usage=usage
        )
        # Effort and thinking are not exposed by `claude -p`, so they must not
        # appear in the key as though they were honoured.
        variant = "probe/claude-cli"
    else:
        complete = anthropic_completer(
            model=args.model,
            effort=args.effort,
            thinking=args.thinking,
            max_tokens=args.max_tokens,
            usage=usage,
        )
        # The scaffold slot in the cache key must capture everything that
        # changes the reply, or a re-run at a different effort replays the old
        # answer and prints a spotless null result.
        variant = f"probe/{args.effort}/{args.thinking}"
    cache = NullCache() if args.no_cache else PromptCache(args.cache)

    def run(probe):
        prompt = probe.prompt()
        cached = cache.get(args.model, variant, prompt)
        if cached is None:
            cached = complete(prompt)
            cache.put(args.model, variant, prompt, cached)
        return score(probe, cached)

    print(
        f"scoring {len(probes)} probes on {args.model} "
        f"(effort={args.effort}, thinking={args.thinking})"
    )
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        results = list(pool.map(run, probes))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8") as stream:
        for result in results:
            stream.write(json.dumps(result, ensure_ascii=False) + "\n")

    print()
    print(f"{'probe':<18} {'n':>5} {'exact':>8} {'jaccard':>9}")
    print("-" * 43)
    for kind in PROBE_KINDS:
        group = [r for r in results if r["kind"] == kind]
        if not group:
            continue
        print(
            f"{kind:<18} {len(group):>5} "
            f"{fmean(r['exact'] for r in group):>7.1%} "
            f"{fmean(r['jaccard'] for r in group):>9.3f}"
        )
    print("-" * 43)
    print(
        f"{'ALL':<18} {len(results):>5} "
        f"{fmean(r['exact'] for r in results):>7.1%} "
        f"{fmean(r['jaccard'] for r in results):>9.3f}"
    )
    print(f"\nusage: {usage.format()}")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
