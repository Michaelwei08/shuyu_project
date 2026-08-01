from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

from .game import ObservedMove, battle_outcome
from .types import Owner, PieceKind, Position

MOVABLE_KINDS = frozenset(kind for kind in PieceKind if kind.movable)


@dataclass
class OpponentKnowledge:
    owner: Owner
    possible: dict[Position, frozenset[PieceKind]] = field(default_factory=dict)
    #: Enemy ranks we can *prove* are dead. Most kills are anonymous, but a few
    #: are not, and those few are the ones worth counting -- mines above all,
    #: because a surviving mine is what makes a rear row impassable.
    destroyed: Counter[PieceKind] = field(default_factory=Counter)

    def observe(self, event: ObservedMove) -> None:
        opponent = self.owner.other
        source, target = event.move.src, event.move.dst
        if event.attacker_owner != opponent and event.had_battle:
            self._count_certain_kill(event)
        if event.attacker_owner == opponent:
            candidates = self.possible.pop(source, MOVABLE_KINDS)
            candidates = candidates & MOVABLE_KINDS
            if not event.had_battle:
                self.possible[target] = candidates
                return
            if event.outcome is not None and event.outcome > 0:
                if event.own_kind is None:
                    raise ValueError("己方守子军衔缺失")
                self.possible[target] = frozenset(
                    kind
                    for kind in candidates
                    if battle_outcome(kind, event.own_kind) > 0
                )
            else:
                self.possible.pop(target, None)
            return

        candidates = self.possible.get(target, frozenset(PieceKind))
        if not event.had_battle:
            return
        if event.outcome is not None and event.outcome < 0:
            if event.own_kind is None:
                raise ValueError("己方攻子军衔缺失")
            self.possible[target] = frozenset(
                kind
                for kind in candidates
                if battle_outcome(event.own_kind, kind) < 0
            )
        else:
            self.possible.pop(target, None)

    def _count_certain_kill(self, event: ObservedMove) -> None:
        """Record the enemy deaths we can actually name.

        Most captures are anonymous -- beating a piece with a colonel says only
        that it was weaker. Two cases are not:

        * an **engineer that wins** can only have beaten a mine, since a mine
          is the single rank it defeats;
        * an **engineer that trades** met a bomb or another engineer, so a bomb
          is at worst over-counted by one, which is the safe direction.

        Everything else is left to the casualty estimate.
        """
        if event.own_kind != PieceKind.ENGINEER or event.outcome is None:
            return
        if event.outcome > 0:
            self.destroyed[PieceKind.MINE] += 1
        elif event.outcome == 0:
            self.destroyed[PieceKind.BOMB] += 1

    def forget_missing(self, occupied: set[Position]) -> None:
        self.possible = {
            position: kinds
            for position, kinds in self.possible.items()
            if position in occupied
        }
