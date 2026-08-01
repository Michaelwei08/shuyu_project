# junqi_training (not part of the Ash science)

This folder is a **generated copy** of the training subset of
`D:\Stanford\research\own\fun\military`, parked here only because this repo is
the pipe to the CPU box. It has nothing to do with viral sequencing -- do not
read it as part of this project, and do not edit it here.

Source of truth is `military/`. Regenerate with:

    python scripts/sync_remote.py push

Pure Python, no third-party dependencies, but it needs **Python 3.10+**
(dataclass slots, `zip(strict=)`).

## First: pin an interpreter

`python3` on a shared box is usually the system one and too old. Do not rely on
`conda activate` -- with no argument it selects base, which on some machines
resolves to `/usr` and is not a conda env at all (`conda list` then fails with
"Not a conda environment: /usr").

    conda env list                # NOT `conda list`
    ls -d ~/.conda/envs/*/

Pick one and pin it for the session, then check it before anything else:

    PY=~/.conda/envs/junqi/bin/python     # substitute the env you found
    $PY -c 'import sys; print(sys.version)'
    cd junqi_training && $PY -m unittest discover -s tests

Expect 35 tests, 1 skipped (the web-weights sync check; `web/` is not shipped
here). Pick a worker count too -- this is CPU-bound, and on a dedicated box
most cores is right:

    W=$(( $(nproc) - 2 ))

## Where the 2026-08-01 round landed (806 paired games each)

    blind pricing       +0.0218 +/- 0.0226, p = 0.17   inconclusive, kept
    defender supply     +0.0008 +/- 0.0097, p = 0.47   dead, now ships at 0
    screen cap 2 vs 3   -0.0985 +/- 0.0217             wrong, reverted to 3

Two of those three were my ideas and the measurement rejected them, which is
the harness working. `eval_hq_supply` is not just unproven -- an SE of 0.0097
excludes anything above +0.016. The screen cap was the largest effect in the
round, in the wrong direction, and instrumenting 208 games showed the flag
falling at the same rate under both caps (23/104 vs 21/104), so the stated
reason for capping was not the mechanism either.

## The one run still worth doing

Only the blind-attack pricing is live, and 806 games cannot settle it: the
minimum detectable improvement at that size is +0.037, larger than the effect
itself. Two independent runs put it at +0.0173 and +0.0218, both positive.
Resolving it needs about three times the sample:

    $PY -m junqi.benchmark --games 2400 --seeds 3 --no-history --workers $W \
        --baseline models/ab/pre-2026-08-01.json

At ~2.6 core-seconds per game that is roughly 3.5 core-hours, so a few minutes
of wall time on 48 cores. `--seeds 3` splits it into disjoint opening blocks,
so a verdict does not rest on one lucky set of deployments.

Accept on `p < 0.05` with a positive paired difference. If it rejects, read the
"minimum detectable improvement" line before concluding: above the plausible
effect size it means the sample was too small, not that the change is dead.

`--no-history` matters. The archived models in `models/history/` predate
`blind_battle`, so they would load with the new default and are not the
opponents they were trained as.

Only after that verdict should training resume from the new baseline:

    $PY -m junqi.training --generations 8 --screen-games 200 \
        --accept-games 600 --workers $W

Only `models/` is meant to travel back. After the run, from `military/`:

    python scripts/sync_remote.py pull
    python -m junqi.web_export

Numbers from before 2026-08-01 are not comparable to numbers after it: the
scoring function changed, both sides of every pool game changed with it, and
`make_opening` now draws each side from its own RNG stream, so every opening
differs. Re-baseline rather than comparing against an old printout.
