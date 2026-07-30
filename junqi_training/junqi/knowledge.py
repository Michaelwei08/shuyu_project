from __future__ import annotations

from dataclasses import dataclass, field

from .game import ObservedMove, battle_outcome
from .types import Owner, PieceKind, Position

MOVABLE_KINDS = frozenset(kind for kind in PieceKind if kind.movable)


@dataclass
class OpponentKnowledge:
    owner: Owner
    possible: dict[Position, frozenset[PieceKind]] = field(default_factory=dict)

    def observe(self, event: ObservedMove) -> None:
        opponent = self.owner.other
        source, target = event.move.src, event.move.dst
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

    def forget_missing(self, occupied: set[Position]) -> None:
        self.possible = {
            position: kinds
            for position, kinds in self.possible.items()
            if position in occupied
        }
