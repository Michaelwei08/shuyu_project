"""The deployment strategy space, and the harness that solves over it.

Kept out of `test_game.py` because these guard one idea: that a deployment
family is a *strategy*, so it has to be legal, has to be selectable per side
without disturbing the other side, and has to be mixable.
"""

from __future__ import annotations

import random
import unittest
from collections import Counter

from junqi.arena import DEFAULT_SEARCH, Job, make_opening, play_match
from junqi.bot import BotWeights
from junqi.deployment import (
    FAMILIES,
    headquarters,
    parse_mixture,
    resolve_family,
    strategic_deployment,
    validate_deployment,
)
from junqi.deployment_game import mirror_spec, regret_matching
from junqi.types import Owner, PieceKind


class FamilyTests(unittest.TestCase):
    def test_every_family_deploys_legally_for_both_sides(self) -> None:
        for name in FAMILIES:
            for owner in Owner:
                for seed in range(40):
                    layout = strategic_deployment(
                        owner, random.Random(seed), family=name
                    )
                    self.assertEqual(
                        validate_deployment(layout, owner),
                        [],
                        f"{name} produced an illegal army for {owner}",
                    )

    def test_each_family_actually_does_what_it_claims(self) -> None:
        """A family that silently falls back to `standard` measures nothing."""
        seen: dict[str, Counter] = {}
        doors: dict[str, Counter] = {}
        for name in (
            "standard",
            "decoy-mine",
            "decoy-bomb",
            "screen2",
            "bomb-home",
            "seal-75",
            "seal-always",
        ):
            decoys: Counter = Counter()
            mines: Counter = Counter()
            for seed in range(120):
                layout = strategic_deployment(
                    Owner.BOT, random.Random(seed), family=name
                )
                squares = sorted(headquarters(Owner.BOT))
                flag = next(
                    square
                    for square in squares
                    if layout[square].kind == PieceKind.FLAG
                )
                decoy = next(square for square in squares if square != flag)
                decoys[layout[decoy].kind] += 1
                mines[
                    sum(
                        1
                        for square, piece in layout.items()
                        if piece.kind == PieceKind.MINE
                        and abs(square[0] - flag[0]) + abs(square[1] - flag[1]) == 1
                    )
                ] += 1
            seen[name] = decoys
            doors[name] = mines

        self.assertEqual(set(seen["decoy-mine"]), {PieceKind.MINE})
        self.assertEqual(set(seen["decoy-bomb"]), {PieceKind.BOMB})
        self.assertEqual(
            set(seen["standard"]), {PieceKind.LIEUTENANT, PieceKind.CAPTAIN}
        )
        # Spending a mine on the decoy leaves only two for the flag's doors.
        self.assertEqual(set(doors["decoy-mine"]), {2})
        # The known-bad family must still be the known-bad family.
        self.assertLessEqual(max(doors["screen2"]), 2)
        # Shipped since 2026-08-04: the seal is unconditional, worth
        # +0.0384 +/- 0.0067 over 2400 paired games. Every opening, not 75%.
        self.assertEqual(set(doors["standard"]), {3})
        self.assertEqual(set(doors["seal-always"]), {3})
        # And the superseded generator must still be selectable, or the result
        # above stops being re-checkable.
        self.assertEqual(set(doors["seal-75"]), {2, 3})

    def test_a_family_changes_only_its_own_army(self) -> None:
        """The property the whole paired design rests on."""
        for owner in Owner:
            for name in FAMILIES:
                base = make_opening(11, {owner: "standard"})
                other = make_opening(11, {owner: name})
                mine = {
                    square: piece
                    for square, piece in base.items()
                    if piece.owner == owner.other
                }
                theirs = {
                    square: piece
                    for square, piece in other.items()
                    if piece.owner == owner.other
                }
                self.assertEqual(mine, theirs, f"{name} perturbed the {owner.other} army")

    def test_a_bare_screen_cap_still_works(self) -> None:
        """`--screen-cap 2` predates families and models/ab/ still speaks it."""
        self.assertEqual(resolve_family(2).screen_cap, 2)
        self.assertEqual(resolve_family(None).name, "standard")
        with self.assertRaises(KeyError):
            resolve_family("no-such-family")

    def test_a_mixture_draws_per_game_and_stays_legal(self) -> None:
        family = parse_mixture("mix:standard=1,decoy-mine=3")
        self.assertAlmostEqual(dict(family.mixture)["decoy-mine"], 0.75)
        decoys: Counter = Counter()
        for seed in range(400):
            layout = strategic_deployment(
                Owner.BOT, random.Random(seed), family="mix:standard=1,decoy-mine=3"
            )
            self.assertEqual(validate_deployment(layout, Owner.BOT), [])
            squares = sorted(headquarters(Owner.BOT))
            flag = next(
                square for square in squares if layout[square].kind == PieceKind.FLAG
            )
            decoy = next(square for square in squares if square != flag)
            decoys[layout[decoy].kind == PieceKind.MINE] += 1
        share = decoys[True] / sum(decoys.values())
        self.assertGreater(share, 0.65)
        self.assertLess(share, 0.85)

    def test_a_mixture_is_reproducible_from_the_seed(self) -> None:
        mix = "mix:standard=0.5,uniform=0.5"
        self.assertEqual(make_opening(5, {Owner.BOT: mix}), make_opening(5, {Owner.BOT: mix}))


class MatrixTests(unittest.TestCase):
    def test_the_two_seats_of_a_mirror_cell_play_the_same_game(self) -> None:
        """Why the matrix diagonal is 0.5 exactly, and only needs checking.

        Both seats hold the same weights and the same search budget, and the
        per-seat RNG seed depends on the seat rather than on who is called the
        subject -- so a diagonal cell is one game scored twice, once from each
        side. If this ever stops holding, `A[j][i] = 1 - A[i][j]` stops holding
        with it and the whole upper-triangle shortcut is wrong.
        """
        model = "models/bot_weights.json"
        weights = BotWeights.load(model)
        spec = mirror_spec(model)
        self.assertEqual(
            (spec.samples, spec.beam_width, spec.reply_width), DEFAULT_SEARCH
        )
        results = [
            play_match(
                Job(
                    weights,
                    spec,
                    901_234,
                    side,
                    subject_deployment="standard",
                    opponent_deployment="standard",
                )
            ).result
            for side in (int(Owner.HUMAN), int(Owner.BOT))
        ]
        self.assertAlmostEqual(sum(results) / 2, 0.5, places=9)

    def test_regret_matching_solves_a_game_with_a_known_answer(self) -> None:
        """Rock-paper-scissors: the only Nash is uniform."""
        names = ["rock", "paper", "scissors"]
        beats = {("rock", "scissors"), ("scissors", "paper"), ("paper", "rock")}
        matrix = {
            (a, b): 0.5 if a == b else (1.0 if (a, b) in beats else 0.0)
            for a in names
            for b in names
        }
        strategy = regret_matching(names, matrix, iterations=20_000)
        for share in strategy:
            self.assertAlmostEqual(share, 1 / 3, places=2)

    def test_regret_matching_finds_a_dominant_pure_strategy(self) -> None:
        """When one family beats everything, the answer must not be a mixture."""
        names = ["good", "bad", "worse"]
        matrix = {
            ("good", "good"): 0.5,
            ("good", "bad"): 0.7,
            ("good", "worse"): 0.8,
            ("bad", "good"): 0.3,
            ("bad", "bad"): 0.5,
            ("bad", "worse"): 0.6,
            ("worse", "good"): 0.2,
            ("worse", "bad"): 0.4,
            ("worse", "worse"): 0.5,
        }
        strategy = regret_matching(names, matrix, iterations=20_000)
        self.assertGreater(strategy[0], 0.99)


if __name__ == "__main__":
    unittest.main()
