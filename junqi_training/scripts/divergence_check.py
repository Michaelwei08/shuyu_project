"""Does a coefficient change any games at all? Run this before believing a null.

A paired null means one of two very different things -- the change was tried and
did nothing, or the change was never really tried -- and `compare()` cannot tell
them apart, because both print a small mean. The tell is the *standard error*:
identical games contribute exactly zero to the paired variance, so a
suspiciously precise comparison is evidence of a change that rarely fires, not
of a well-measured one.

    python scripts/divergence_check.py --weight eval_hq_entomb --scales 0.5 1 2

Worked example, 2026-08-03. `eval_hq_entomb` measured +0.0029 +/- 0.0032 over
2400 games -- an SE a fifth of the 0.0142 the deployment families produced on
the same harness. That looked like excellent power. It was dilution:

    scale   diverged   mean diff        SD
      0.5       8.4%     +0.0032    0.1130
      1.0      11.2%     +0.0082    0.1540
      2.0      26.9%     -0.0003    0.2296

Two things only this table shows. The term fires on a tenth of games at the
scale that was tested, so the aggregate is nine parts untouched game. And
raising it to 2.0 fires on 2.4x as many games and nets *exactly* nothing, which
is what kills the dilution hypothesis -- if the term were good on the games it
touched, touching more of them would show it. Louder is not better here; it
starts suppressing the probes that were correct, since a headquarters attack
takes the flag 60% of the time.

Do not read the conditional effect (mean / diverged) as the adoption criterion.
It is a diagnostic: divergence is caused by the treatment, so conditioning on it
is not a post-hoc slice, but the aggregate is still what D012 names as the bar.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import fields
from pathlib import Path
from statistics import fmean, stdev

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from junqi.arena import build_jobs, run_jobs, seeds_for  # noqa: E402
from junqi.bot import BotWeights  # noqa: E402
from junqi.opponents import standard_pool  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=Path("models/bot_weights.json"))
    parser.add_argument("--anchor", type=Path, default=Path("models/defaults.json"))
    parser.add_argument("--weight", required=True, help="要检验的系数名")
    parser.add_argument("--scales", type=float, nargs="+", default=[0.5, 1.0, 2.0])
    parser.add_argument("--games", type=int, default=320, help="每档的对局数")
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--offset", type=int, default=880_000)
    arguments = parser.parse_args()

    known = {descriptor.name for descriptor in fields(BotWeights)}
    if arguments.weight not in known:
        raise SystemExit(f"no such coefficient: {arguments.weight}")

    shipped = BotWeights.load(arguments.model)
    pool = standard_pool(anchor=str(arguments.anchor) if arguments.anchor else None)
    seeds = seeds_for(arguments.games, pool, offset=arguments.offset)

    baseline = build_jobs(shipped, pool, seeds)
    arms = []
    for scale in arguments.scales:
        loud = BotWeights.load(arguments.model)
        setattr(loud, arguments.weight, scale)
        arms.append((scale, build_jobs(loud, pool, seeds)))

    jobs = baseline + [job for _, arm in arms for job in arm]
    print(
        f"{arguments.weight}: baseline {getattr(shipped, arguments.weight)} "
        f"vs {arguments.scales}",
    )
    print(f"{len(jobs)} games ({len(baseline)} per arm, {len(arms) + 1} arms)",
          flush=True)
    played = run_jobs(jobs, arguments.workers)
    width = len(baseline)
    base = played[:width]

    print(
        f"\n{'scale':>7} {'games':>7} {'diverged':>10} {'mean diff':>11} "
        f"{'SD':>8} {'SE':>8} {'per touched':>12}"
    )
    print("-" * 68)
    for index, (scale, _) in enumerate(arms, start=1):
        window = played[index * width : (index + 1) * width]
        differences = [
            candidate.score - incumbent.score
            for candidate, incumbent in zip(window, base, strict=True)
            if candidate is not None and incumbent is not None
        ]
        if not differences:
            print(f"{scale:>7} every paired match failed")
            continue
        moved = [value for value in differences if abs(value) > 1e-12]
        share = len(moved) / len(differences)
        deviation = stdev(differences) if len(differences) > 1 else 0.0
        mean = fmean(differences)
        print(
            f"{scale:>7} {len(differences):>7} {share:>9.1%} {mean:>+11.4f} "
            f"{deviation:>8.4f} "
            f"{deviation / len(differences) ** 0.5:>8.4f} "
            f"{(mean / share if share else float('nan')):>+12.4f}"
        )
    print(
        "\nRead the divergence column first. A coefficient that fires on a tenth "
        "of games cannot be judged by an aggregate that is nine parts untouched "
        "game -- and if raising the scale fires on more games without moving the "
        "mean, the term is neutral rather than under-powered."
    )


if __name__ == "__main__":
    main()
