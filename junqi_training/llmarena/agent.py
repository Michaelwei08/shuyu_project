"""An LLM as a drop-in ``choose_move`` player.

The interface is the engine's existing duck type -- ``choose_move(game, owner)``
returning a legal ``Move`` -- so an ``LLMAgent`` slots into ``arena.play_match``
next to ``RandomBot`` and ``SearchBot`` with no special-casing.

**The raw proposal and the repair are recorded separately.** ``choose_move``
always returns something legal, and the transcript keeps what the model actually
asked for. That is what lets one set of paid calls score more than one penalty
regime: a *forfeit* rule is re-derived offline by truncating each game at its
first illegal proposal (see :func:`first_illegal_ply`), so it costs nothing
extra. *Retry* is the exception and cannot be derived -- an extra attempt draws
a fresh sample and sends the game down a different path -- so it has to be run
as its own condition.
"""

from __future__ import annotations

import random
import re
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field
from typing import Any

from junqi.game import Game
from junqi.types import Move, Owner, parse_position

from .belief import BeliefTracker
from .cache import NullCache, PromptCache
from .view import LEGAL_BLOCK_MARKER, Scaffold, build_observation, render

Completer = Callable[[str], str]

#: Tolerant on the way in, strict on the way out. Accepts ``A10-A9``,
#: ``A10 -> A9``, ``a10 a9``; the *last* match in the response wins, so a model
#: may reason about other moves first and still answer cleanly on its last line.
MOVE_PATTERN = re.compile(
    r"([A-Ea-e])\s*(1[0-2]|[1-9])\s*[\s\-—–>→]+\s*([A-Ea-e])\s*(1[0-2]|[1-9])"
)


def parse_response(text: str) -> Move | None:
    """Pull the intended move out of a completion, or ``None``."""
    matches = MOVE_PATTERN.findall(text or "")
    if not matches:
        return None
    src_col, src_row, dst_col, dst_row = matches[-1]
    try:
        return Move(
            parse_position(f"{src_col}{src_row}"),
            parse_position(f"{dst_col}{dst_row}"),
        )
    except ValueError:
        return None


@dataclass
class TurnRecord:
    """One decision, enough of it to re-score offline."""

    ply: int
    side: str
    scaffold: str
    legal_count: int
    #: Every raw completion, in order. Length > 1 only under the retry policy.
    responses: list[str] = field(default_factory=list)
    #: What the model asked for on its *first* attempt, or None if unparseable.
    proposed: str | None = None
    #: Whether that first proposal was legal.
    proposal_legal: bool = False
    #: The move actually applied to the board.
    played: str = ""
    #: How the played move was arrived at: "model" | "retry" | "random".
    source: str = "model"

    def to_dict(self) -> dict[str, Any]:
        return {
            "ply": self.ply,
            "side": self.side,
            "scaffold": self.scaffold,
            "legal_count": self.legal_count,
            "responses": self.responses,
            "proposed": self.proposed,
            "proposal_legal": self.proposal_legal,
            "played": self.played,
            "source": self.source,
        }


def first_illegal_ply(transcript: Sequence[TurnRecord]) -> int | None:
    """Ply at which a forfeit rule would have ended the game, if ever.

    The whole point of separating proposal from repair: this reads a run that
    was played under ``fallback`` and answers what ``forfeit`` would have done,
    without another call.
    """
    for record in transcript:
        if not record.proposal_legal:
            return record.ply
    return None


class LLMAgent:
    """Plays by prompting a completion function.

    ``complete`` is injected rather than constructed here so the whole class is
    testable with a stub and no network.
    """

    def __init__(
        self,
        complete: Completer,
        *,
        model: str,
        scaffold: Scaffold,
        cache: PromptCache | NullCache | None = None,
        repair: str = "fallback",
        max_retries: int = 2,
        seed: int | None = None,
        transcript: list[TurnRecord] | None = None,
        cache_namespace: str | None = None,
    ) -> None:
        if repair not in ("fallback", "retry"):
            raise ValueError(
                f"unknown repair policy {repair!r}; "
                "forfeit is derived offline, not played live"
            )
        self.complete = complete
        self.model = model
        self.scaffold = scaffold
        self.cache = cache if cache is not None else NullCache()
        #: Cache key namespace. Must capture *everything* that changes the
        #: reply, not just the scaffold: effort and thinking do too, and a run
        #: at a new effort that replayed the old answers would look like a
        #: perfectly clean null result.
        self.cache_namespace = (
            scaffold.name if cache_namespace is None else cache_namespace
        )
        self.repair = repair
        self.max_retries = max_retries
        self.rng = random.Random(seed)
        self.transcript: list[TurnRecord] = (
            transcript if transcript is not None else []
        )
        self.tracker: BeliefTracker | None = None

    def _ask(self, prompt: str, variant: int) -> str:
        cached = self.cache.get(self.model, self.cache_namespace, prompt, variant)
        if cached is not None:
            return cached
        response = self.complete(prompt)
        self.cache.put(self.model, self.cache_namespace, prompt, response, variant)
        return response

    def choose_move(self, game: Game, owner: Owner | None = None) -> Move:
        player = game.turn if owner is None else owner
        legal = game.legal_moves(player)
        if not legal:
            raise ValueError("当前玩家没有合法走法")

        if self.tracker is None or self.tracker.owner != player:
            self.tracker = BeliefTracker(player)
        belief = self.tracker.update(game)

        observation = build_observation(game, player, belief)
        prompt = render(observation, self.scaffold)
        legal_set = set(legal)

        record = TurnRecord(
            ply=game.move_count,
            side=player.name,
            scaffold=self.scaffold.name,
            legal_count=len(legal),
        )
        self.transcript.append(record)

        attempts = 1 if self.repair == "fallback" else 1 + self.max_retries
        for attempt in range(attempts):
            text = self._ask(
                prompt if attempt == 0 else self._nudge(prompt), attempt
            )
            record.responses.append(text)
            move = parse_response(text)
            if attempt == 0:
                record.proposed = None if move is None else str(move)
                record.proposal_legal = move is not None and move in legal_set
            if move is not None and move in legal_set:
                record.played = str(move)
                record.source = "model" if attempt == 0 else "retry"
                return move

        fallback = self.rng.choice(legal)
        record.played = str(fallback)
        record.source = "random"
        return fallback

    @staticmethod
    def _nudge(prompt: str) -> str:
        return (
            prompt
            + "\n\n上一次的回答不是一个合法走法。请从上面列出的走法中重新选择，"
            "并确保最后一行只有 起点-终点 这一种格式。"
        )


# --- completers -----------------------------------------------------------
# Day 1 ships only stubs. A real provider is a function of the same shape, so
# adding one touches this registry and nothing else.


def scripted_completer(responses: Sequence[str]) -> Completer:
    """Replays a fixed list, then repeats the last entry."""
    queue = list(responses)

    def complete(_prompt: str) -> str:
        if not queue:
            return responses[-1] if responses else ""
        return queue.pop(0)

    return complete


def random_legal_completer(seed: int | None = None) -> Completer:
    """Answers with a legal move lifted out of the prompt itself.

    Only usable with a scaffold that lists legal moves, which is the point: it
    exercises the full prompt/parse/apply path at zero cost and gives the rest
    of the harness a control player whose scaffolding sensitivity is known to
    be exactly nil.
    """
    rng = random.Random(seed)

    def complete(prompt: str) -> str:
        # Read only the legal-move block. The move log and the format example
        # in the closing instruction also look like moves, and picking one of
        # those would make this control player illegal at random -- which is
        # the one thing it must not be.
        start = prompt.find(LEGAL_BLOCK_MARKER)
        block = (
            prompt[start:].split("\n\n", 1)[0] if start >= 0 else prompt
        )
        moves = MOVE_PATTERN.findall(block)
        if not moves:
            return "无法作答"
        col_a, row_a, col_b, row_b = rng.choice(moves)
        return f"我选择 {col_a}{row_a}-{col_b}{row_b}"

    return complete


COMPLETERS: dict[str, Callable[[dict[str, str]], Completer]] = {
    "random-legal": lambda options: random_legal_completer(
        seed=int(options.get("seed", "0"))
    ),
}


def make_completer(name: str, options: dict[str, str]) -> Completer:
    if name == "claude-cli" and name not in COMPLETERS:
        from .cli_completer import register as register_cli

        register_cli()
    if name == "anthropic" and name not in COMPLETERS:
        # Registered on demand so importing this module never pulls in the SDK,
        # and so a worker process that only unpickles an AgentSpec still
        # resolves the backend without the caller remembering to register it.
        from .anthropic_completer import register

        register()
    if name not in COMPLETERS:
        raise ValueError(
            f"unknown completer {name!r}; available: {sorted(COMPLETERS)}"
        )
    return COMPLETERS[name](options)


def build_agent(spec: Any, _weights: Any, seed: int) -> LLMAgent:
    """Factory behind ``AgentSpec(builder="llmarena.agent:build_agent")``.

    Runs inside the worker process, so the completion client is constructed
    there rather than pickled across.
    """
    from .view import SCAFFOLDS

    options = dict(spec.options)
    scaffold_name = options.get("scaffold", "legal")
    if scaffold_name not in SCAFFOLDS:
        raise ValueError(f"unknown scaffold {scaffold_name!r}")
    cache_root = options.get("cache")
    effort = options.get("effort", "-")
    thinking = options.get("thinking", "-")
    return LLMAgent(
        make_completer(options.get("completer", "random-legal"), options),
        model=options.get("model", "stub"),
        scaffold=SCAFFOLDS[scaffold_name],
        cache=PromptCache(cache_root) if cache_root else NullCache(),
        cache_namespace=f"{scaffold_name}/{effort}/{thinking}",
        repair=options.get("repair", "fallback"),
        max_retries=int(options.get("max_retries", "2")),
        seed=seed,
    )
