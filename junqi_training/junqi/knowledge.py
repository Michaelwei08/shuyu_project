from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

from .game import ObservedMove, battle_outcome
from .types import PIECE_COUNTS, Owner, PieceKind, Position

MOVABLE_KINDS = frozenset(kind for kind in PieceKind if kind.movable)
ALL_KINDS = frozenset(PieceKind)


def _dead_attacker_candidates(
    defender: PieceKind, outcome: int
) -> frozenset[PieceKind]:
    """Ranks an enemy *attacker* could have been, given that it died.

    ``outcome`` is signed from the attacker's view, so the attacker died when it
    is negative (we held) or zero (mutual). The attacker moved, so it is
    movable -- never a mine, never the flag.
    """
    return frozenset(
        kind
        for kind in MOVABLE_KINDS
        if (battle_outcome(kind, defender) < 0 and outcome < 0)
        or (battle_outcome(kind, defender) == 0 and outcome == 0)
    )


def _dead_defender_candidates(
    attacker: PieceKind, outcome: int
) -> frozenset[PieceKind]:
    """Ranks an enemy *defender* could have been, given that it died.

    The defender died when the attacker won (>0) or both fell (0).

    **The flag stays in.** Excluding it looks safe -- capturing the flag ends
    the game, so surely it never matters -- and it is the bug this deduction
    shipped with for one commit. The final record of a won game *is* a flag
    capture, and `observe` sees it like any other. Since an engineer beats
    exactly the flag and mines, dropping the flag turned that record into the
    singleton `{MINE}`, and three engineer wins then "proved" mines extinct
    while one was still on the board. Caught by
    `test_an_extinct_rank_really_has_no_survivors` on seed 14.

    Keeping the flag costs the engineer-win mine pin (now `{FLAG, MINE}`, not a
    singleton) and keeps the chain that matters: a mine can only trade with a
    bomb, and that is an attacker-side deduction the flag cannot reach.
    """
    return frozenset(
        kind
        for kind in ALL_KINDS
        if (battle_outcome(attacker, kind) > 0 and outcome > 0)
        or (battle_outcome(attacker, kind) == 0 and outcome == 0)
    )


@dataclass
class OpponentKnowledge:
    owner: Owner
    possible: dict[Position, frozenset[PieceKind]] = field(default_factory=dict)
    #: Enemy ranks we can *prove* are dead. Most kills are anonymous, but a few
    #: are not, and those few are the ones worth counting -- mines above all,
    #: because a surviving mine is what makes a rear row impassable.
    destroyed: Counter[PieceKind] = field(default_factory=Counter)
    #: One entry per enemy piece known to have died, holding the ranks it could
    #: have been. Feeds :meth:`eliminate_dead_ranks`; costs nothing until that
    #: is called, so tracking it is unconditional and using it is gated.
    dead_enemy: list[frozenset[PieceKind]] = field(default_factory=list)
    #: Read an enemy railway corner as an engineer. Only engineers turn corners,
    #: so this is a certainty rather than a prior -- and it is the only rank
    #: information a *quiet* move carries, which is why this deduction is the one
    #: thing `observe` used to discard. Gated by `use_engineer_deduction`,
    #: because the pool opponents run this class too and a plain `if` would apply
    #: to both sides of a paired comparison and cancel.
    deduce_engineers: bool = False

    def observe(self, event: ObservedMove) -> None:
        opponent = self.owner.other
        source, target = event.move.src, event.move.dst
        if event.had_battle:
            self._record_enemy_death(event)
        if event.attacker_owner != opponent and event.had_battle:
            self._count_certain_kill(event)
        if event.attacker_owner == opponent:
            candidates = self.possible.pop(source, MOVABLE_KINDS)
            candidates = candidates & MOVABLE_KINDS
            if self.deduce_engineers and event.engineer_only:
                # Certain, so it overrides rather than narrows. An empty
                # intersection would mean an earlier deduction was wrong, and the
                # topology is the thing we are surest of.
                narrowed = candidates & frozenset({PieceKind.ENGINEER})
                candidates = narrowed or frozenset({PieceKind.ENGINEER})
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

    def _record_enemy_death(self, event: ObservedMove) -> None:
        """Note the rank-range of any enemy piece that just died."""
        if event.outcome is None or event.own_kind is None:
            return
        if event.attacker_owner == self.owner.other:
            # They attacked us. Their attacker died on a loss or a trade.
            if event.outcome <= 0:
                candidates = _dead_attacker_candidates(event.own_kind, event.outcome)
            else:
                return
        else:
            # We attacked them. Their defender died on our win or a trade.
            if event.outcome >= 0:
                candidates = _dead_defender_candidates(event.own_kind, event.outcome)
            else:
                return
        if candidates:
            self.dead_enemy.append(candidates)

    def eliminate_dead_ranks(
        self, known_alive: frozenset[PieceKind] = frozenset()
    ) -> frozenset[PieceKind]:
        """Remove ranks the casualty record proves are extinct.

        ``observe`` reasons one square at a time: it can say "whatever survived
        there beats a colonel" but never "that rank is already all dead, so it
        cannot be here". This closes that gap by counting, and it is the
        deduction a language model made unprompted on this position set -- both
        bombs accounted for, therefore the piece that traded with our general
        *was* the general, therefore no survivor is one.

        The rule is a fixpoint over two sound steps:

        1. A death whose candidate set is a singleton names a rank exactly. The
           common source is our own mine: nothing but a bomb can trade with a
           mine, so a mine that goes down mutually kills a bomb, provably.
        2. Once as many deaths are pinned to rank R as the army contains, R is
           extinct -- so R leaves every live square's belief *and* every other
           casualty's candidate set, which can pin those exactly, and round
           again.

        Sound but not complete: it does not solve the full bipartite matching
        between casualties and the roster, so it under-reports rather than
        over-reports. Returns the extinct ranks.
        """
        # A rank we can prove is still standing cannot be a casualty. The one
        # such proof this game offers is free and load-bearing: the flag is
        # revealed exactly when the commander dies, so an unrevealed flag means
        # the commander lives -- which turns every `{COMMANDER, BOMB}` trade
        # into a pinned bomb. That is the step the language model used and the
        # engine was missing; without it the deduction never fires at all
        # (measured: 0 extinctions in 60 games).
        pools = [
            set(candidates) - known_alive if len(candidates - known_alive) else set(candidates)
            for candidates in self.dead_enemy
        ]
        extinct: set[PieceKind] = set()
        while True:
            pinned: Counter[PieceKind] = Counter()
            for pool in pools:
                if len(pool) == 1:
                    pinned[next(iter(pool))] += 1
            fresh = {
                kind
                for kind, count in pinned.items()
                if count >= PIECE_COUNTS[kind] and kind not in extinct
            }
            if not fresh:
                break
            extinct |= fresh
            for pool in pools:
                if len(pool) > 1:
                    # Safe only because `fresh` ranks are fully accounted for by
                    # *other* casualties; never empty a pool on a bad inference.
                    narrowed = pool - fresh
                    if narrowed:
                        pool.clear()
                        pool |= narrowed
        if extinct:
            self.possible = {
                position: kinds - extinct for position, kinds in self.possible.items()
            }
        return frozenset(extinct)

    def forget_missing(self, occupied: set[Position]) -> None:
        self.possible = {
            position: kinds
            for position, kinds in self.possible.items()
            if position in occupied
        }
