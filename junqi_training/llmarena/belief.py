"""Incremental upkeep of what one side has legally deduced.

This is a re-implementation of ``SearchBot._update_knowledge``, not a call into
it: that method is bound to a ``SearchBot`` (which also builds a ``HeuristicBot``
and a search RNG) and an LLM agent needs the belief without the search. The
duplication is guarded by ``test_the_tracker_agrees_with_the_search_bot``, which
plays a game and asserts the two ``possible`` maps stay identical ply by ply --
so a change to the engine's deduction fails here rather than silently giving the
LLM a different picture than the bot it is being compared against.

Every read is leak-free: piece owners, the ``revealed`` bit, and the public move
log. Never a hidden rank.
"""

from __future__ import annotations

from junqi.game import Game
from junqi.knowledge import OpponentKnowledge
from junqi.types import Owner, PieceKind, Position


class BeliefTracker:
    """Folds the public move log into a per-square set of possible ranks."""

    def __init__(self, owner: Owner, eliminate_dead_ranks: bool = False) -> None:
        self.owner = owner
        self.knowledge = OpponentKnowledge(owner)
        self.processed_records = 0
        #: Mirrors `BotWeights.use_rank_elimination`. Default off so the probe
        #: labels keep matching what the shipped bot actually deduces -- turning
        #: it on here without turning it on in the engine would silently make
        #: the battery measure a bot that does not exist.
        self.eliminate_dead_ranks = eliminate_dead_ranks

    def update(self, game: Game) -> dict[Position, frozenset[PieceKind]]:
        """Consume any new records and return the current belief."""
        for event in game.observations(self.owner, self.processed_records):
            self.knowledge.observe(event)
        self.processed_records = len(game.records)

        opponent_positions = {
            position
            for position, piece in game.board.items()
            if piece.owner == self.owner.other
        }
        self.knowledge.forget_missing(opponent_positions)

        if self.eliminate_dead_ranks:
            self.knowledge.eliminate_dead_ranks()

        if commander_dead(game, self.owner.other):
            # Their commander is gone, so no surviving piece can be one.
            self.knowledge.possible = {
                position: kinds - {PieceKind.COMMANDER}
                for position, kinds in self.knowledge.possible.items()
            }
        return self.knowledge.possible


def commander_dead(game: Game, side: Owner) -> bool:
    """Whether ``side`` has lost its commander.

    A commander's death is exactly what reveals that side's flag, so the
    ``revealed`` bit on a headquarters square is a public proxy for it. Reads
    occupancy and ``revealed`` only.
    """
    return any(game.board[square].revealed for square in game.flag_candidates(side))
