from __future__ import annotations

import random
from collections import Counter

from .board import HEADQUARTERS
from .bot import BotWeights, HeuristicBot, _distance, _piece_value
from .deployment import headquarters, rear_rows
from .game import Game
from .knowledge import OpponentKnowledge
from .types import PIECE_COUNTS, Move, Owner, Piece, PieceKind, Position


class SearchBot:
    """Belief-sampling search that never scores an unknown piece by its true rank."""

    def __init__(
        self,
        weights: BotWeights | None = None,
        seed: int | None = None,
        samples: int = 6,
        beam_width: int = 10,
        reply_width: int = 4,
    ) -> None:
        self.weights = weights or BotWeights()
        self.rng = random.Random(seed)
        self.base_seed = 0 if seed is None else seed
        self.samples = samples
        self.beam_width = beam_width
        self.reply_width = reply_width
        self.heuristic = HeuristicBot(self.weights, seed=seed)
        self.knowledge: OpponentKnowledge | None = None
        self.processed_records = 0

    def choose_move(self, game: Game, owner: Owner | None = None) -> Move:
        player = game.turn if owner is None else owner
        self._update_knowledge(game, player)
        self.heuristic.knowledge = (
            None if self.knowledge is None else self.knowledge.possible
        )
        moves = game.legal_moves(player)
        if not moves:
            raise ValueError("当前玩家没有合法走法")
        ranked = sorted(
            ((self.heuristic._score(game, move, player), move) for move in moves),
            reverse=True,
            key=lambda item: item[0],
        )[: self.beam_width]
        # Common random numbers. Deriving the sampled worlds from
        # (seed, ply) rather than from a running RNG means two candidate weight
        # sets playing the same colour in the same game face byte-identical
        # hidden states, so a comparison measures policy and not luck.
        stream = random.Random(self.base_seed * 1_000_003 + game.move_count)
        sample_seeds = [stream.randrange(2**32) for _ in range(self.samples)]
        scored: list[tuple[float, Move]] = []
        for base_score, move in ranked:
            rollout = sum(
                self._rollout(game, move, player, seed) for seed in sample_seeds
            ) / len(sample_seeds)
            scored.append((base_score + rollout, move))
        best = max(score for score, _ in scored)
        return self.rng.choice(
            [move for score, move in scored if score >= best - 1e-9]
        )

    def _update_knowledge(self, game: Game, owner: Owner) -> None:
        if self.knowledge is None or self.knowledge.owner != owner:
            self.knowledge = OpponentKnowledge(owner)
            self.processed_records = 0
        for event in game.observations(owner, self.processed_records):
            self.knowledge.observe(event)
        self.processed_records = len(game.records)
        opponent_positions = {
            position
            for position, piece in game.board.items()
            if piece.owner == owner.other
        }
        self.knowledge.forget_missing(opponent_positions)

    def _rollout(
        self, game: Game, move: Move, owner: Owner, sample_seed: int
    ) -> float:
        sampled = self._determinize(game, owner.other, random.Random(sample_seed))
        sampled.turn = owner
        sampled.apply(move)
        if sampled.winner == owner:
            return self.weights.eval_terminal
        if sampled.winner == owner.other:
            return -self.weights.eval_terminal
        opponent = owner.other
        replies = sampled.legal_moves(opponent)
        if not replies:
            return self._state_value(sampled, owner)

        # Take the worst case over the opponent's most plausible answers rather
        # than trusting a single greedy reply.
        reply_bot = HeuristicBot(self.weights, seed=sample_seed ^ 0xA5A5A5A5)
        ranked = sorted(
            (
                (reply_bot._score(sampled, reply, opponent, quick=True), reply)
                for reply in replies
            ),
            reverse=True,
            key=lambda item: item[0],
        )[: self.reply_width]
        worst = None
        for _, reply in ranked:
            child = sampled.clone()
            child.turn = opponent
            child.apply(reply)
            value = self._state_value(child, owner)
            if worst is None or value < worst:
                worst = value
        return worst if worst is not None else self._state_value(sampled, owner)

    def _determinize(
        self, game: Game, hidden_owner: Owner, rng: random.Random
    ) -> Game:
        sampled = game.clone()
        hidden_positions = [
            position
            for position, piece in sampled.board.items()
            if piece.owner == hidden_owner and not piece.revealed
        ]
        revealed = Counter(
            piece.kind
            for piece in sampled.board.values()
            if piece.owner == hidden_owner and piece.revealed
        )
        for position in hidden_positions:
            sampled.board.pop(position)

        constraints = self.knowledge.possible if self.knowledge is not None else {}
        assignment: dict[Position, PieceKind] | None = None
        for _ in range(8):
            hidden_kinds = self._sample_survivors(
                len(hidden_positions), revealed, rng
            )
            assignment = self._assign_constrained(
                hidden_owner, hidden_positions, hidden_kinds, constraints, rng
            )
            if assignment is not None:
                break
        if assignment is None:
            # Belief and sampled survivors cannot be reconciled; fall back to an
            # unconstrained shuffle so the search still has a world to play in.
            shuffled = self._sample_survivors(len(hidden_positions), revealed, rng)
            rng.shuffle(shuffled)
            assignment = dict(zip(hidden_positions, shuffled, strict=True))
        for position, kind in assignment.items():
            sampled.board[position] = Piece(hidden_owner, kind)
        return sampled

    @staticmethod
    def _sample_survivors(
        count: int, revealed: Counter[PieceKind], rng: random.Random
    ) -> list[PieceKind]:
        """Guess which ranks are still alive without looking at the real board.

        Battles are anonymous, so the bot does not know which of the opponent's
        ranks died -- it may only subtract the ranks it has actually seen.
        """
        pool = Counter(PIECE_COUNTS)
        pool.subtract(revealed)
        survivors = list(pool.elements())
        for _ in range(max(0, len(survivors) - count)):
            # The flag is alive while the game is running, so kill it last --
            # but never return more kinds than there are squares to fill.
            removable = [
                index
                for index, kind in enumerate(survivors)
                if kind != PieceKind.FLAG
            ] or list(range(len(survivors)))
            survivors.pop(rng.choice(removable))
        return survivors

    def _assign_constrained(
        self,
        owner: Owner,
        positions: list[Position],
        kinds: list[PieceKind],
        constraints: dict[Position, frozenset[PieceKind]],
        rng: random.Random,
    ) -> dict[Position, PieceKind] | None:
        remaining = Counter(kinds)
        ordered = positions.copy()
        rng.shuffle(ordered)
        all_kinds = frozenset(PieceKind)
        ordered.sort(key=lambda position: len(constraints.get(position, all_kinds)))
        assignment: dict[Position, PieceKind] = {}

        def allowed(position: Position, kind: PieceKind) -> bool:
            if kind == PieceKind.FLAG and position not in headquarters(owner):
                return False
            if kind == PieceKind.MINE and position[0] not in rear_rows(owner):
                return False
            return kind in constraints.get(position, all_kinds)

        budget = [len(ordered) * 40]

        def assign(index: int) -> bool:
            if index == len(ordered):
                return True
            budget[0] -= 1
            if budget[0] < 0:
                return False
            position = ordered[index]
            pending = ordered[index + 1 :]
            choices = [
                kind
                for kind, count in remaining.items()
                if count > 0 and allowed(position, kind)
            ]
            rng.shuffle(choices)
            # Try the ranks with the fewest legal homes left first. Flags and
            # mines only fit a handful of squares, and filling those squares
            # with anything else first is what makes this search thrash.
            choices.sort(
                key=lambda kind: sum(1 for slot in pending if allowed(slot, kind))
            )
            for kind in choices:
                remaining[kind] -= 1
                assignment[position] = kind
                if assign(index + 1):
                    return True
                remaining[kind] += 1
                assignment.pop(position)
            return False

        return assignment if assign(0) else None

    def _state_value(self, game: Game, owner: Owner) -> float:
        weights = self.weights
        if game.winner == owner:
            return weights.eval_terminal
        if game.winner == owner.other:
            return -weights.eval_terminal
        material = sum(
            _piece_value(piece.kind) * (1 if piece.owner == owner else -1)
            for piece in game.board.values()
        )
        mobility = len(game.legal_moves(owner)) - len(game.legal_moves(owner.other))
        return (
            material * weights.eval_material
            + mobility * weights.eval_mobility
            + self._flag_pressure(game, owner)
        )

    def _flag_pressure(self, game: Game, owner: Owner) -> float:
        """Reward closing on the enemy headquarters, punish losing our own.

        Without this the evaluation is pure material and the bot has no reason
        to ever go for the win.
        """
        weights = self.weights
        value = 0.0
        attacking = game.flag_candidates(owner.other)
        reach = self._closest_raider(game, owner, attacking)
        if reach is not None:
            value += (12 - reach) * (
                weights.eval_hq_attack_certain
                if len(attacking) == 1
                else weights.eval_hq_attack
            )
        defending = game.flag_candidates(owner)
        threat = self._closest_raider(game, owner.other, defending)
        if threat is not None:
            value -= (12 - threat) * (
                weights.eval_hq_defense_certain
                if len(defending) == 1
                else weights.eval_hq_defense
            )
            if threat <= 1:
                # An enemy standing next to a headquarters that may hold our
                # flag takes it next ply unless we remove it. The linear term
                # above barely distinguishes that from a distant piece.
                value -= weights.eval_hq_breach
        value += weights.eval_hq_guard * self._guards(game, owner, defending)
        return value

    @staticmethod
    def _guards(game: Game, owner: Owner, squares: list[Position]) -> int:
        """Own movable pieces standing next to a headquarters we must hold.

        Without this, stepping off the square that shields the flag cost only
        the `protect_flag` distance term -- 0.4 against a ~24 point capture.
        """
        if not squares:
            return 0
        return sum(
            1
            for position, piece in game.board.items()
            if piece.owner == owner
            and piece.kind.movable
            and position not in HEADQUARTERS
            and any(_distance(position, square) == 1 for square in squares)
        )

    @staticmethod
    def _closest_raider(
        game: Game, side: Owner, targets: list[Position]
    ) -> int | None:
        if not targets:
            return None
        return min(
            (
                _distance(position, target)
                for position, piece in game.board.items()
                if piece.owner == side
                and piece.kind.movable
                and position not in HEADQUARTERS
                for target in targets
            ),
            default=None,
        )
