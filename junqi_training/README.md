# junqi_training (not part of the Ash science)

This folder is a **generated copy** of the training subset of
`D:\Stanford\research\own\fun\military`, parked here only because this repo is
the pipe to the CPU box. It has nothing to do with viral sequencing -- do not
read it as part of this project, and do not edit it here.

Source of truth is `military/`. Regenerate with:

    python scripts/sync_remote.py push

On the compute box:

    cd junqi_training
    python3 -m junqi.training --generations 8 --screen-games 200 \
        --accept-games 600 --workers <vCPUs>

Only `models/` is meant to travel back. After the run, from `military/`:

    python scripts/sync_remote.py pull
    python -m junqi.web_export

Pure Python, no third-party dependencies. Python 3.10+.
