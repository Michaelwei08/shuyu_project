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

## Acceptance runs for the 2026-08-01 changes

Four changes need a paired verdict. Each is isolated by a baseline model in
`models/ab/`, because `compare()` varies weights and not code. At ~2.6
core-seconds per game, 800 games is roughly 35 core-minutes per run.

    # 1. blind-attack pricing + defender supply, together
    $PY -m junqi.benchmark --games 800 --no-history --workers $W \
        --baseline models/ab/pre-2026-08-01.json

    # 2. blind-attack pricing alone
    $PY -m junqi.benchmark --games 800 --no-history --workers $W \
        --baseline models/ab/no-blind.json

    # 3. defender supply alone
    $PY -m junqi.benchmark --games 800 --no-history --workers $W \
        --baseline models/ab/no-supply.json

Deployment is not a weight, so it has its own switch. `--screen-cap` changes
only the subject's own layout -- the opposing army on a given seed is
untouched, which is what keeps the comparison paired:

    # 4. two-mine flag screen vs the old three-mine seal
    $PY -m junqi.benchmark --games 800 --no-history --workers $W \
        --baseline models/bot_weights.json \
        --screen-cap 2 --baseline-screen-cap 3

Accept on `p < 0.05` with a positive paired difference. If a run rejects, check
the printed "minimum detectable improvement" before concluding anything: above
the plausible effect size it means the sample was too small, not that the
change is dead. Add `--seeds 3` for a release gate.

`--no-history` matters here. The archived models in `models/history/` predate
`blind_battle` and `eval_hq_supply`, so they would be loaded with the new
defaults and are not the opponents they were trained as.

Only after those verdicts should training resume from the new baseline:

    $PY -m junqi.training --generations 8 --screen-games 200 \
        --accept-games 600 --workers $W

Only `models/` is meant to travel back. After the run, from `military/`:

    python scripts/sync_remote.py pull
    python -m junqi.web_export

Numbers from before 2026-08-01 are not comparable to numbers after it: the
scoring function changed, both sides of every pool game changed with it, and
`make_opening` now draws each side from its own RNG stream, so every opening
differs. Re-baseline rather than comparing against an old printout.
