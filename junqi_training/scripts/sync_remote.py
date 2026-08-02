"""Ship the training subset to a git repo that the compute box pulls from.

`military/` is not a git repo, and `web/` inside it is a separate one with its
own remote, so the whole folder is awkward to move (and 691 MB of it is
`node_modules`). Training needs about a quarter of a megabyte.

    python scripts/sync_remote.py push          # military -> transport repo
    python scripts/sync_remote.py pull          # trained models -> military

`military/` stays the single source of truth for code: `push` overwrites the
copy in the transport repo, never the other way round. Only `models/` travels
back, because that is the only thing the remote produces.
"""

from __future__ import annotations

import argparse
import filecmp
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REPO = Path(r"D:\Stanford\research\Ash")
FOLDER = "junqi_training"

# `scripts/` travels too: `generalisation_check.py` is now a required step
# before adopting anything fitted, so the compute box needs it.
CODE = ["junqi", "tests", "scripts"]
FILES = ["pyproject.toml"]
MODELS = "models"
SKIP = {"__pycache__", ".pytest_cache"}

README = """\
# junqi_training (not part of the Ash science)

This folder is a **generated copy** of the training subset of
`D:\\Stanford\\research\\own\\fun\\military`, parked here only because this repo is
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

Expect 42 tests, 1 skipped (the web-weights sync check; `web/` is not shipped
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

## Blind pricing: closed, dead

Five measurements, converging on zero:

    n=78,    corrupt anchor           +0.19            pure noise in hindsight
    n=806 x2,corrupt anchor           +0.0173, +0.0218
    n=2418,  corrupt anchor           +0.0066 +/- 0.0133, p=0.31
    n=806 x3,repaired, scales 2/4/6   -0.0132, -0.0074, +0.0077, all p>0.35
    n=2418,  repaired anchor          -0.0014 +/- 0.0118, p=0.55

Repairing the anchor did NOT revive it -- notable, because the same repair
doubled the value-function result. The instrument was hiding the leaf
evaluation; it was not hiding this.

It still ships at blind_battle = 9.0 on the same grounds use_move_distance
ships at 1 despite measuring +0.0021: the old term is not weaker, it is
BACKWARDS -- it scored the commander lowest and the engineer highest for the
identical blind attack, against observed outcomes of 3W/1T/0L and 1W/2T/5L.
Equal against this pool; only one is correct. That is a judgement, not a result.

Do not spend more games on it. models/ab/blind-{2,4,6}.json are kept only so the
negative can be reproduced.

## Two instrument repairs -- read before quoting any old number

**The anchor was not a yardstick.** models/defaults.json was missing seven
fields, and BotWeights.load back-fills missing keys with TODAY's dataclass
defaults -- so the anchor moved on every commit that added a coefficient, with
the file untouched. It had drifted into running unknown_risk = 0.12 AND
blind_battle = 9.0 at once: both blind-attack pricings, used by ten of thirteen
pool opponents, making the pool over-cautious about attacking unknowns. That is
the exact axis a human beat the bot on (76% vs 46%). Repaired; all 34 fields are
now explicit and a test fails if a new one is added without rerunning:

    python scripts/rebuild_anchor.py --write

**Paired SEs are optimistic by ~12%.** compare() treats matches as independent
but 26 share each opening. Measured design effect 1.26 over 806 games. Multiply
any reported SE by 1.12. It changes none of the 2026-08-01 conclusions.

## The learned leaf evaluation -- ADOPTED

junqi/value.py fits the leaf evaluation to self-play outcomes: 96 parameters,
stdlib SGD, weighted by eval_value_scale, now shipping at 15.0.

    500 games,  n=400,  corrupt anchor    -0.0013 +/- 0.0318   nothing
    2000 games, n=400,  corrupt anchor    +0.0631 +/- 0.0303   screen only
    2000 games, n=1600, corrupt anchor    +0.0210 +/- 0.0153   not significant
    1600/3200/6400, n=2400, corrupt       +0.0293/+0.0373/+0.0288  all p<0.01
    1600 games, n=2400, REPAIRED anchor   +0.0654 +/- 0.0123, p ~ 1e-6

All checks are against four opponents excluded from the model's training. Two
things to carry forward: an n=400 screen means nothing here, and the effect
roughly DOUBLED once the anchor was repaired -- the instrument was hiding it.

The data curve is flat above ~1600 games (the three sizes differ by under one
SE), and capacity is not the constraint either (6->30 features gains +0.038
held-out AUC, 30->96 gains nothing, train AUC below test AUC throughout). So
neither more games nor a bigger model is the next lever.

To re-baseline or to challenge the adoption:

    $PY -m junqi.benchmark --games 2400 --seeds 3 --no-history --workers $W \
        --baseline models/ab/value-off.json

To refit after any change to the anchor or the weights:

    $PY -m junqi.value_training --games 1600 --workers $W --stride 4 \
        --exclude search-deep search-mid selective selective-strict
    $PY scripts/generalisation_check.py --workers $W --games 2400 --scales 15

Only after that verdict should training resume from the new baseline:

    $PY -m junqi.training --generations 8 --screen-games 200 \\
        --accept-games 600 --workers $W

Only `models/` is meant to travel back. After the run, from `military/`:

    python scripts/sync_remote.py pull
    python -m junqi.web_export

Numbers from before 2026-08-01 are not comparable to numbers after it: the
scoring function changed, both sides of every pool game changed with it, and
`make_opening` now draws each side from its own RNG stream, so every opening
differs. Re-baseline rather than comparing against an old printout.
"""


def _copy_tree(source: Path, destination: Path) -> int:
    copied = 0
    if destination.exists():
        shutil.rmtree(destination)
    for item in source.rglob("*"):
        if any(part in SKIP for part in item.parts):
            continue
        target = destination / item.relative_to(source)
        if item.is_dir():
            target.mkdir(parents=True, exist_ok=True)
        else:
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(item, target)
            copied += 1
    return copied


def push(repo: Path, with_models: bool = False) -> None:
    destination = repo / FOLDER
    destination.mkdir(parents=True, exist_ok=True)
    total = 0
    for name in CODE:
        total += _copy_tree(ROOT / name, destination / name)
    for name in FILES:
        shutil.copy2(ROOT / name, destination / name)
        total += 1

    existing = destination / MODELS
    if with_models or not existing.exists():
        # `_copy_tree` deletes the destination first, so doing this on every
        # code sync would silently destroy weights the compute box trained but
        # has not pushed back yet. Models travel remote -> local by default.
        total += _copy_tree(ROOT / MODELS, existing)
        print("models: overwritten from local")
    else:
        print("models: left alone (pass --with-models to overwrite)")

    (destination / "README.md").write_text(README, encoding="utf-8", newline="\n")
    print(f"pushed {total + 1} files -> {destination}")
    print(f"  cd {repo} && git add {FOLDER} && git commit -m 'Sync junqi training subset'")


def pull(repo: Path) -> None:
    source = repo / FOLDER / MODELS
    if not source.exists():
        raise SystemExit(f"no models in transport repo: {source}")
    destination = ROOT / MODELS
    changed: list[str] = []
    for item in sorted(source.rglob("*.json")):
        target = destination / item.relative_to(source)
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists() and filecmp.cmp(item, target, shallow=False):
            continue
        shutil.copy2(item, target)
        changed.append(str(target.relative_to(ROOT)))
    if not changed:
        print("models already up to date")
        return
    print("updated:")
    for name in changed:
        print(f"  {name}")
    print("\nnow run: python -m junqi.web_export   (otherwise the web bot is stale)")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("push", "pull"))
    parser.add_argument("--repo", type=Path, default=DEFAULT_REPO)
    parser.add_argument(
        "--with-models",
        action="store_true",
        help="也用本地 models/ 覆盖传输仓库（会丢弃远程尚未回传的训练结果）",
    )
    arguments = parser.parse_args()
    if not (arguments.repo / ".git").exists():
        raise SystemExit(f"not a git repo: {arguments.repo}")
    if arguments.action == "push":
        push(arguments.repo, arguments.with_models)
    else:
        pull(arguments.repo)


if __name__ == "__main__":
    main()
