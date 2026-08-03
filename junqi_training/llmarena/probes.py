"""A dense diagnostic battery with exact ground truth.

Full matches cost ~42 calls and return one bit. These probes cost one call and
return a graded answer against a label the engine computes exactly, which is the
difference between a weekend of measurement and a month of it. They also split
the loss into layers that a win rate cannot:

* ``legal_moves`` -- rule execution. Rail sliding, engineer turning, the river
  gate on columns A/C/E. Label: ``Game._destinations``.
* ``flag_candidates`` -- public deduction. Which headquarters can still hold the
  enemy flag, from occupancy alone. Label: ``Game.flag_candidates``. This is the
  one inference the whole game hinges on and it needs no hidden information at
  all, so failing it is a clean, strong result.
* ``belief`` -- private inference. Which ranks a square can still hold given the
  anonymous battle results so far. Label: ``OpponentKnowledge.possible``.

**Positions are stratified, not sampled uniformly.** Measured over ten games,
only ~1.4 squares per position carry a non-trivial deduced rank set, so uniform
sampling produces a battery whose labels are mostly empty. Generation is cheap
(~2.6 core-seconds per game), so it over-generates and keeps only the positions
where the label discriminates.

Each probe carries a serialised :class:`Observation` rather than a
``(seed, ply)`` pointer, so the battery is a self-contained artifact: it does
not drift when the bot weights that generated it change.
"""

from __future__ import annotations

import json
import random
import re
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from junqi.bot import BotWeights, HeuristicBot
from junqi.game import Game
from junqi.types import SYMBOLS, Owner, Position, format_position

from .belief import BeliefTracker
from .view import MOVABLE_KINDS_COUNT, Observation, Scaffold, build_observation, render

PROBE_KINDS = ("legal_moves", "flag_candidates", "belief")

#: Each probe is rendered with the scaffold flag that would hand over its own
#: answer switched **off**. Enforced by
#: ``test_a_probe_prompt_never_contains_its_own_answer``.
PROBE_SCAFFOLDS: dict[str, Scaffold] = {
    "legal_moves": Scaffold(
        "probe-legal_moves", legal_moves=False, flag_candidates=True, belief=True
    ),
    "flag_candidates": Scaffold(
        "probe-flag_candidates", legal_moves=True, flag_candidates=False, belief=True
    ),
    "belief": Scaffold(
        "probe-belief", legal_moves=True, flag_candidates=True, belief=False
    ),
}

_SYMBOL_TO_KIND = {symbol: kind for kind, symbol in SYMBOLS.items()}
_SQUARE_PATTERN = re.compile(r"\b([A-Ea-e])\s*(1[0-2]|[1-9])\b")


@dataclass(frozen=True)
class Probe:
    """One question with an exact answer."""

    probe_id: str
    kind: str
    observation: Observation
    #: Free-text question appended after the rendered position.
    question: str
    #: Canonical answer, as a sorted list of strings (squares or rank names).
    label: tuple[str, ...]
    meta: dict[str, Any]

    def prompt(self) -> str:
        return render(self.observation, PROBE_SCAFFOLDS[self.kind], self.question)

    def to_dict(self) -> dict[str, Any]:
        return {
            "probe_id": self.probe_id,
            "kind": self.kind,
            "observation": self.observation.to_dict(),
            "question": self.question,
            "label": list(self.label),
            "meta": self.meta,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "Probe":
        return cls(
            probe_id=data["probe_id"],
            kind=data["kind"],
            observation=Observation.from_dict(data["observation"]),
            question=data["question"],
            label=tuple(data["label"]),
            meta=data.get("meta", {}),
        )


# --- answer parsing -------------------------------------------------------


def parse_squares(text: str) -> set[str]:
    """Coordinates named in the final answer line."""
    line = _answer_line(text)
    return {
        f"{column.upper()}{row}" for column, row in _SQUARE_PATTERN.findall(line)
    }


def parse_kinds(text: str) -> set[str]:
    """Rank names in the final answer line, as ``PieceKind`` names."""
    line = _answer_line(text)
    return {
        _SYMBOL_TO_KIND[symbol].name
        for symbol in _SYMBOL_TO_KIND
        if symbol in line
    }


def _answer_line(text: str) -> str:
    lines = [line.strip() for line in (text or "").splitlines() if line.strip()]
    return lines[-1] if lines else ""


def score(probe: Probe, response: str) -> dict[str, Any]:
    """Grade one answer. Exact set match is primary, Jaccard is the gradient."""
    predicted = (
        parse_kinds(response)
        if probe.kind == "belief"
        else parse_squares(response)
    )
    truth = set(probe.label)
    union = predicted | truth
    return {
        "probe_id": probe.probe_id,
        "kind": probe.kind,
        "exact": predicted == truth,
        "jaccard": len(predicted & truth) / len(union) if union else 1.0,
        "predicted": sorted(predicted),
        "label": sorted(truth),
    }


# --- generation -----------------------------------------------------------


def _legal_probe(
    obs: Observation, rng: random.Random, probe_id: str
) -> Probe | None:
    """Ask for one piece's destinations, preferring a piece with rail options."""
    by_source: dict[Position, list[Position]] = {}
    for move in obs.legal:
        by_source.setdefault(move.src, []).append(move.dst)
    # A piece with one or two road steps tests nothing; rail sliding is where
    # the rules actually bite.
    interesting = [
        square for square, targets in by_source.items() if len(targets) >= 5
    ]
    if not interesting:
        return None
    square = rng.choice(sorted(interesting))
    kind = next(piece.kind for piece in obs.own if piece.square == square)
    targets = tuple(sorted(format_position(t) for t in by_source[square]))
    return Probe(
        probe_id=probe_id,
        kind="legal_moves",
        observation=obs,
        question=(
            f"问题：你在 {format_position(square)} 的{SYMBOLS[kind]}，"
            "这一步可以走到哪些格子？请列出全部落点，不要遗漏也不要多写。\n"
            "最后一行只输出坐标，用空格分隔，例如：A5 B6 C7"
        ),
        label=targets,
        meta={"square": format_position(square), "kind": kind.name},
    )


def _flag_probe(obs: Observation, probe_id: str) -> Probe | None:
    if not obs.enemy_flag_squares:
        return None
    return Probe(
        probe_id=probe_id,
        kind="flag_candidates",
        observation=obs,
        question=(
            "问题：对方的军旗现在可能在哪些格子？请列出全部可能位置。\n"
            "最后一行只输出坐标，用空格分隔，例如：B1 D1"
        ),
        label=tuple(sorted(format_position(s) for s in obs.enemy_flag_squares)),
        meta={"count": len(obs.enemy_flag_squares)},
    )


def _belief_probe(
    obs: Observation, rng: random.Random, probe_id: str
) -> Probe | None:
    informative = [
        (square, kinds)
        for square, kinds in obs.belief
        if 0 < len(kinds) < MOVABLE_KINDS_COUNT
    ]
    if not informative:
        return None
    square, kinds = rng.choice(informative)
    return Probe(
        probe_id=probe_id,
        kind="belief",
        observation=obs,
        question=(
            f"问题：根据上面的行棋记录，对方在 {format_position(square)} 的那枚棋子，"
            "军衔还有哪些可能？请列出全部仍然可能的军衔。\n"
            "最后一行只输出军衔名称，用空格分隔，例如：司 军 师"
        ),
        label=tuple(sorted(kind.name for kind in kinds)),
        meta={"square": format_position(square), "size": len(kinds)},
    )


def generate(
    seeds: Sequence[int],
    weights: BotWeights | None = None,
    max_plies: int = 120,
    sample_every: int = 6,
) -> Iterator[Probe]:
    """Play throwaway games and harvest positions where a label discriminates.

    Uses ``HeuristicBot`` on both sides: the battery is about positions, not
    about strength, and the heuristic is roughly twenty times cheaper than the
    search bot per ply.
    """
    from junqi.arena import make_opening

    weights = weights or BotWeights()
    for seed in seeds:
        game = Game(board=make_opening(seed), turn=Owner(seed % 2))
        players = {
            owner: HeuristicBot(weights, seed=seed * 2 + int(owner))
            for owner in Owner
        }
        trackers = {owner: BeliefTracker(owner) for owner in Owner}
        rng = random.Random(seed * 7_919 + 11)

        while not game.over and game.move_count < max_plies:
            mover = game.turn
            belief = trackers[mover].update(game)
            if game.move_count % sample_every == 0:
                obs = build_observation(game, mover, belief)
                stem = f"s{seed}-p{game.move_count}-{mover.name.lower()}"
                for probe in (
                    _legal_probe(obs, rng, f"legal_moves-{stem}"),
                    _flag_probe(obs, f"flag_candidates-{stem}"),
                    _belief_probe(obs, rng, f"belief-{stem}"),
                ):
                    if probe is not None:
                        yield probe
            game.apply(players[mover].choose_move(game))


def balanced(
    probes: Sequence[Probe], per_kind: int, seed: int = 0
) -> list[Probe]:
    """Take an equal number of each kind, spread across games.

    ``flag_candidates`` is additionally balanced between the two-candidate case
    and the decided one-candidate case, which is the half that actually matters
    and is much rarer.
    """
    rng = random.Random(seed)
    chosen: list[Probe] = []
    for kind in PROBE_KINDS:
        pool = [probe for probe in probes if probe.kind == kind]
        if kind == "flag_candidates":
            narrow = [p for p in pool if p.meta.get("count", 2) == 1]
            wide = [p for p in pool if p.meta.get("count", 2) != 1]
            rng.shuffle(narrow)
            rng.shuffle(wide)
            half = per_kind // 2
            picked = narrow[:half] + wide[: per_kind - len(narrow[:half])]
        else:
            picked = list(pool)
            rng.shuffle(picked)
            picked = picked[:per_kind]
        chosen.extend(picked)
    return chosen


def write_jsonl(probes: Sequence[Probe], path: str | Path) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8") as stream:
        for probe in probes:
            stream.write(json.dumps(probe.to_dict(), ensure_ascii=False) + "\n")
    return destination


def read_jsonl(path: str | Path) -> list[Probe]:
    with Path(path).open("r", encoding="utf-8") as stream:
        return [Probe.from_dict(json.loads(line)) for line in stream if line.strip()]
