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

Blind pricing at `blind_battle = 9.0` then failed too, at 2418 paired games:
**+0.0066 +/- 0.0133, p = 0.31**, minimum detectable +0.0218. The point estimate
shrank as the sample grew (+0.19 at 78 games, +0.017 and +0.022 at 806, +0.0066
at 2418), which is noise regressing to zero.

## The one thing still open

That run split cleanly by opponent class rather than randomly:

    helped  defensive +4.5  hqrush +3.2  heuristic +3.0  engineer +2.7
            material +2.7   hqrush-careful +2.2  heuristic-quiet +1.6
            selective +0.9  random +0.5
    hurt    search-mid -5.4  search-deep -4.6  search-shallow -3.0
            selective-strict -1.6

All eight greedy opponents positive, all three search opponents negative. The
likely cause is that the fix changed two things at once: the rank *ordering*
(which really was inverted) and the *magnitude*. The branch spans 14.6 points
across the rank order at 9.0 where the old term spanned 0.73, so it drowns out
the rollout the search is supposed to be steering by.

`models/ab/blind-{2,4,6}.json` keep the ordering and turn the volume down
(spans 3.2 / 6.5 / 9.7). Screen all three -- about 130s each:

    for s in 2 4 6; do
      $PY -m junqi.benchmark --games 800 --no-history --workers $W \
          --model models/ab/blind-$s.json \
          --baseline models/ab/pre-2026-08-01.json
    done

These are screens and picking the best of three inflates the false-positive
rate, so a winner has to be confirmed on a larger sample before it means
anything:

    $PY -m junqi.benchmark --games 2400 --seeds 3 --no-history --workers $W \
        --model models/ab/blind-<best>.json \
        --baseline models/ab/pre-2026-08-01.json

If none of them clears `p < 0.05` at 2400 games, the whole line is dead: set
`blind_battle` to 0, restore `unknown_risk`, and the answer is that the pool
cannot distinguish these pricings at all.

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
