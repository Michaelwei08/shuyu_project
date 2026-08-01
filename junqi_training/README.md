# junqi_training (not part of the Ash science)

This folder is a **generated copy** of the training subset of
`D:\Stanford\research\own\fun\military`, parked here only because this repo is
the pipe to the CPU box. It has nothing to do with viral sequencing -- do not
read it as part of this project, and do not edit it here.

Source of truth is `military/`. Regenerate with:

    python scripts/sync_remote.py push

Pure Python, no third-party dependencies. Python 3.10+. Start with:

    cd junqi_training
    python3 -m unittest discover -s tests

## Acceptance runs for the 2026-08-01 changes

Three changes need a paired verdict. Each is isolated by a baseline model in
`models/ab/`, because `compare()` varies weights and not code. Set `--workers`
to your vCPU count; at ~2.6 core-seconds per game, 800 games is roughly 35
core-minutes per run.

    # 1. blind-attack pricing + defender supply, together
    python3 -m junqi.benchmark --games 800 --no-history --workers N \
        --baseline models/ab/pre-2026-08-01.json

    # 2. blind-attack pricing alone
    python3 -m junqi.benchmark --games 800 --no-history --workers N \
        --baseline models/ab/no-blind.json

    # 3. defender supply alone
    python3 -m junqi.benchmark --games 800 --no-history --workers N \
        --baseline models/ab/no-supply.json

Deployment is not a weight, so it has its own switch. `--screen-cap` changes
only the subject's own layout -- the opposing army on a given seed is
untouched, which is what keeps the comparison paired:

    # 4. two-mine flag screen vs the old three-mine seal
    python3 -m junqi.benchmark --games 800 --no-history --workers N \
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

    python3 -m junqi.training --generations 8 --screen-games 200 \
        --accept-games 600 --workers N

Only `models/` is meant to travel back. After the run, from `military/`:

    python scripts/sync_remote.py pull
    python -m junqi.web_export

Numbers from before 2026-08-01 are not comparable to numbers after it: the
scoring function changed, both sides of every pool game changed with it, and
`make_opening` now draws each side from its own RNG stream, so every opening
differs. Re-baseline rather than comparing against an old printout.
