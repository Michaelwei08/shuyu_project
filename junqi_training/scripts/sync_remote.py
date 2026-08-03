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
# before adopting anything fitted, so the compute box needs it. `llmarena/`
# travels because `tests/` imports it -- shipping the tests without the
# package would break the remote test run, which is the first thing anyone
# does after pulling. Missing entries are skipped rather than fatal, so this
# list can name a package that does not exist yet.
CODE = ["junqi", "llmarena", "tests", "scripts"]
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

## Then: repin the anchor, because `models/` does not travel

`push` deliberately never overwrites `models/`, since this box is what *produces*
trained weights. But `models/defaults.json` is not produced here -- it is a
code-derived yardstick, and `BotWeights.load` fills missing keys from the
dataclass, so **every coefficient added upstream since the last models commit
resolves silently to today's default** and `test_the_anchor_pins_every_
coefficient` fails. That has happened twice (`use_rank_elimination`, then
`eval_hq_entomb`). Both times the repair reported *no behaviour change* -- the
field was already resolving to the value it then pinned -- so this is hygiene,
not a re-baseline. Run it after every pull, before the tests:

    $PY scripts/rebuild_anchor.py --write     # read its "BEHAVIOUR CHANGES" line
    $PY scripts/make_ablations.py             # ablations track *this* box's model
    $PY -m unittest discover -s tests

If `rebuild_anchor` ever prints an actual behaviour change, stop: numbers
measured against the old anchor are no longer comparable and need re-baselining.

One test skips (the web-weights sync check; `web/` is not shipped here);
everything else must pass. Pick a worker count too -- this is CPU-bound, and on a dedicated box
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

## The learned leaf evaluation -- NOT adopted, ships at 0

junqi/value.py fits the leaf evaluation to self-play outcomes: 96 parameters,
stdlib SGD, weighted by eval_value_scale. How it got there:

    500 games,  n=400,  corrupt anchor    -0.0013 +/- 0.0318   nothing
    2000 games, n=400,  corrupt anchor    +0.0631 +/- 0.0303   screen only
    2000 games, n=1600, corrupt anchor    +0.0210 +/- 0.0153   not significant
    1600/3200/6400, n=2400, corrupt       +0.0293/+0.0373/+0.0288  all p<0.01
    1600 games, n=2400, REPAIRED anchor   +0.0654 +/- 0.0123, p ~ 1e-6

SETTLED: real, and too small to adopt. One model (the 6400-game fit),
seed-clustered p in brackets:

    4 held-out opponents, n=2400   +0.0389 +/- 0.0126   p = 0.0010 (0.0029)
    full pool,            n=2418   +0.0161 +/- 0.0108   p = 0.067  (0.092)
    full pool,            n=4836   +0.0093 +/- 0.0078   p = 0.11   (0.14)

Clears the HARDER test -- the one built to catch overfitting to the pool -- and
does not clear the aggregate that D012 names as the criterion. At n=4836 the
95% CI is [-0.008, +0.026] and the minimum detectable improvement is +0.0128,
above the observed +0.0093.

Not contradictory, arithmetically consistent: the held-out four are 4/13 of the
pool, so +0.0389 on them contributes +0.0120 to the aggregate, and the observed
+0.0093 says the other nine average zero. The term helps against the strongest
opponents (heuristic +6.4, search-shallow +5.1) and does nothing against the
rest.

eval_value_scale therefore SHIPS AT 0; models/ab/value-on.json is the
switched-on end. That is a bar-and-cost call, not a verdict of worthless:
+0.8 points of aggregate win rate is not worth 19% of a 420ms browser budget.
Flip it on if the bar changes, if the browser stops binding, or if the pool is
rebuilt so its ceiling is not the subject's own policy class (search-mid).

    $PY -m junqi.benchmark --games 4800 --seeds 3 --no-history --workers $W \
        --model models/ab/value-on.json

TRAP, and it cost a 910-second run: after adoption bot_weights.json CARRIES the
scale under test, so any check comparing a candidate against "the shipped model"
races it against itself and prints +0.0000 +/- 0.0000, p = 1.0000 -- which reads
exactly like a clean negative. Both scripts now zero the scale to build their
baseline and abort on zero variance. An n=400 screen also still means nothing
here; the effect roughly doubled once the anchor was repaired.

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

## The pool now has a ceiling, and it cheats

`search-mid` topped out near 60% against the subject and IS the subject's own
policy class with anchor weights, so the pool could not grade anything past
parity-with-self. Two attempts to fix that; only one worked.

WORKED -- `OracleSearchBot`, same weights but it reads every hidden rank.
Subject wins only 40.6% against it over 240 games, so it is a real ceiling
rather than another near-copy. It is also the CHEAPEST member: at gamma 1 every
sampled world is the same world, so samples collapses to 1. `oracle-half`
(gamma 0.5) is Suphx's dropout schedule and sits between the two. This does not
breach the hidden-rank invariant, which binds the SUBJECT -- the agent that
ships; a test asserts no shipped module references it.

FAILED -- the best-response exploiter. `scripts/exploiter.py` aims the ordinary
training loop at a single opponent, which is what a best-response oracle is.
Against search-mid: **0 of 10 generations accepted, +0.0%** (62.0% before,
62.0% after, 600 games each way). That is worth more than it looks. Coefficient
tuning was already known not to improve the AGGREGATE; this says the weight
space cannot even be bent to beat ONE FIXED OPPONENT better. The 34-dimensional
linear policy class is saturated, and "train harder" is closed as an avenue.

Two bugs that failure exposed, both fixed, both silent:
* `train` saves the unchanged incumbent when nothing is accepted, so the banked
  exploiter was a BYTE-IDENTICAL COPY OF THE SUBJECT -- and discovery would have
  added the bot to its own judging pool. exploiter.py now deletes rather than
  banks a model that did not move or that misses --min-gain.
* `train` derived its anchor from `output.parent`, so writing into
  models/exploiters/ forked a SECOND defaults.json there and trained against it.
  The anchor is now passed explicitly.

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
        source = ROOT / name
        if not source.exists():
            print(f"skipped {name}/ (not in this checkout)")
            continue
        total += _copy_tree(source, destination / name)
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
