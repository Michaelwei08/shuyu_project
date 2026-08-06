from __future__ import annotations

import random
from collections import Counter

from .board import COLUMNS, HEADQUARTERS, ROWS, move_distance
from .bot import BotWeights, HeuristicBot, _distance, _piece_value
from .deployment import headquarters, rear_rows
from .game import Game
from .knowledge import OpponentKnowledge
from .types import PIECE_COUNTS, Move, Owner, Piece, PieceKind, Position
from .value import load_default as load_value_model

# Nothing on this board is more than five moves from anything else, so the old
# `12 - distance` horizon left the term almost flat.
EVAL_HORIZON = 6

#: How close an enemy must get before defender supply starts being charged for.
#: Two moves, because a raider one move out is already too late to reinforce
#: against -- that is exactly the position every replayed loss ended in.
REINFORCE_RANGE = 2
#: Answers we want available. One is not enough: every loss came from spending
#: the only defender and having nothing behind it when the square was refilled.
MIN_COVER = 2


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
            scored.append(
                (base_score * self.weights.search_base_weight + rollout, move)
            )
        best = max(score for score, _ in scored)
        return self.rng.choice(
            [move for score, move in scored if score >= best - 1e-9]
        )

    def _update_knowledge(self, game: Game, owner: Owner) -> None:
        if self.knowledge is None or self.knowledge.owner != owner:
            self.knowledge = OpponentKnowledge(owner)
            self.processed_records = 0
        # Set every call, not just on construction: a weight can differ between
        # two candidates sharing one agent class, and the flag has to follow it.
        self.knowledge.deduce_engineers = bool(self.weights.use_engineer_deduction)
        for event in game.observations(owner, self.processed_records):
            self.knowledge.observe(event)
        self.processed_records = len(game.records)
        opponent_positions = {
            position
            for position, piece in game.board.items()
            if piece.owner == owner.other
        }
        self.knowledge.forget_missing(opponent_positions)
        if self.weights.use_rank_elimination:
            self.knowledge.eliminate_dead_ranks()
        if self._commander_dead(game, owner.other):
            # Their commander is gone, so no surviving piece can be one.
            self.knowledge.possible = {
                position: kinds - {PieceKind.COMMANDER}
                for position, kinds in self.knowledge.possible.items()
            }

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
        reply_bot.knowledge = self._reply_belief(sampled, owner, sample_seed)
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

    def _reply_belief(
        self, sampled: Game, owner: Owner, sample_seed: int
    ) -> dict[Position, frozenset[PieceKind]] | None:
        """What the modelled opponent is assumed to know about *our* army.

        Reads only ``owner``'s own pieces -- our ranks, which we obviously may
        look at -- so this adds no information channel from the hidden side.
        `reply_insight` is the share of our pieces the replier sees; at 0 this
        returns ``None`` and the reply ranking is byte-identical to before.

        Drawn from ``sample_seed`` rather than a running RNG so that two
        candidate weight sets still meet identical worlds and the paired
        comparison stays exact.
        """
        insight = self.weights.reply_insight
        if insight <= 0.0:
            return None
        rng = random.Random(sample_seed ^ 0x5F5F5F5F)
        return {
            position: frozenset({piece.kind})
            for position, piece in sampled.board.items()
            if piece.owner == owner and rng.random() < insight
        }

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
        commander_dead = self._commander_dead(sampled, hidden_owner)
        rear = rear_rows(hidden_owner)
        rear_slots = sum(1 for position in hidden_positions if position[0] in rear)
        assignment: dict[Position, PieceKind] | None = None
        for _ in range(8):
            hidden_kinds = self._sample_survivors(
                len(hidden_positions),
                revealed,
                rng,
                commander_dead,
                self.knowledge.destroyed if self.knowledge is not None else None,
                rear_slots,
            )
            assignment = self._assign_constrained(
                hidden_owner, hidden_positions, hidden_kinds, constraints, rng
            )
            if assignment is not None:
                break
        if assignment is None:
            # Belief and sampled survivors cannot be reconciled. Fall back to a
            # world that ignores the *belief* but still obeys the *rules* -- the
            # old fallback was a bare shuffle, which happily put mines in
            # midfield and the flag outside a headquarters, i.e. worlds the game
            # could never have produced.
            shuffled = self._sample_survivors(
                len(hidden_positions), revealed, rng, rear_slots=rear_slots
            )
            assignment = self._place_by_rules(
                hidden_owner, hidden_positions, shuffled, rng
            )
        for position, kind in assignment.items():
            sampled.board[position] = Piece(hidden_owner, kind)
        return sampled

    @staticmethod
    def _commander_dead(game: Game, side: Owner) -> bool:
        """A commander's death is what reveals its own flag, so the two are
        equivalent. Free, exact, and it collapses a lot of guessing: once the
        enemy commander is gone, whatever beat our major general can only be
        the general."""
        return any(
            game.board[square].revealed for square in game.flag_candidates(side)
        )

    @staticmethod
    def _sample_survivors(
        count: int,
        revealed: Counter[PieceKind],
        rng: random.Random,
        commander_dead: bool = False,
        destroyed: Counter[PieceKind] | None = None,
        rear_slots: int | None = None,
    ) -> list[PieceKind]:
        """Guess which ranks are still alive without looking at the real board.

        Battles are anonymous, so most casualties have to be estimated. But the
        estimate was picking victims uniformly, which quietly killed mines the
        bot had no reason to believe were dead -- and a mine is exactly the
        piece whose survival makes a rear row impassable. So: subtract the
        ranks actually seen, subtract the deaths we can prove, then draw the
        remaining casualties from the ranks that plausibly die.

        With no proven mine kills this leaves all three alive, which is what
        turns "four pieces left in their back rows" into "three mines and the
        flag" once the constraint solver places them.
        """
        pool = Counter(PIECE_COUNTS)
        pool.subtract(revealed)
        if commander_dead:
            pool[PieceKind.COMMANDER] = 0
        for kind, dead in (destroyed or Counter()).items():
            pool[kind] = max(0, pool[kind] - dead)
        if rear_slots is not None:
            # Mines never move, so a live mine must be standing on one of the
            # enemy's own rear-row squares -- there cannot be more of them than
            # there are such squares still occupied. Pure occupancy, no hidden
            # rank read.
            #
            # Without this the estimate keeps all three mines alive forever,
            # because a mine only dies to an engineer (proven, and subtracted
            # above) or to a bomb -- and a bomb trades with *every* rank, so a
            # bomb kill is unprovable and never subtracted. The phantom mines
            # are not merely a bad prior: `_assign_constrained` may then have
            # more mines than legal rear squares to put them on, fail all eight
            # attempts, and fall through to the unconstrained shuffle that
            # scatters mines across the middle of the board.
            # The flag is also stuck in the rear (a headquarters sits there), so
            # it consumes one of the same slots.
            room = max(0, rear_slots - pool[PieceKind.FLAG])
            pool[PieceKind.MINE] = min(pool[PieceKind.MINE], room)

        survivors = list(pool.elements())
        # A flag never dies while the game runs, and a mine is far more likely
        # to be alive than a mobile piece, so both are the last to be guessed
        # dead -- but the cap above has already bounded how many can be here.
        protected = {PieceKind.FLAG, PieceKind.MINE}
        for _ in range(max(0, len(survivors) - count)):
            removable = [
                index
                for index, kind in enumerate(survivors)
                if kind not in protected
            ] or [
                index
                for index, kind in enumerate(survivors)
                if kind != PieceKind.FLAG
            ] or list(range(len(survivors)))
            survivors.pop(rng.choice(removable))
        return survivors

    @staticmethod
    def _place_by_rules(
        owner: Owner,
        positions: list[Position],
        kinds: list[PieceKind],
        rng: random.Random,
    ) -> dict[Position, PieceKind]:
        """Seat every rank on a square the rules actually allow.

        Hardest first: a flag fits only a headquarters and a mine only the rear
        two rows, so those are seated before the ranks that fit anywhere. Used
        only when the belief-constrained solver gives up -- this ignores belief,
        but an illegal world is worse than an uninformed one.
        """
        rear = rear_rows(owner)
        homes = headquarters(owner)
        free = sorted(positions)
        rng.shuffle(free)

        def take(fits) -> Position | None:
            for index, position in enumerate(free):
                if fits(position):
                    return free.pop(index)
            return None

        def difficulty(kind: PieceKind) -> int:
            if kind == PieceKind.FLAG:
                return 0
            return 1 if kind == PieceKind.MINE else 2

        assignment: dict[Position, PieceKind] = {}
        for kind in sorted(kinds, key=difficulty):
            if kind == PieceKind.FLAG:
                position = take(lambda square: square in homes)
            elif kind == PieceKind.MINE:
                position = take(lambda square: square[0] in rear)
            else:
                position = None
            if position is None:
                position = take(lambda _: True)
            if position is None:
                break
            assignment[position] = kind
        return assignment

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
        own_moves = game.legal_moves(owner)
        mobility = len(own_moves) - len(game.legal_moves(owner.other))
        concealment = self._commander_shield(game, owner) - self._commander_shield(
            game, owner.other
        )
        squeeze = weights.eval_immobilize * (
            2.0 ** -self._mobile_count(game, owner.other)
            - 2.0 ** -self._mobile_count(game, owner)
        )
        total = (
            material * weights.eval_material
            + mobility * weights.eval_mobility
            + concealment * weights.eval_commander
            + squeeze
            + self._flag_pressure(game, owner, own_moves)
        )
        if weights.eval_value_scale:
            # The learned leaf evaluation, added rather than substituted so the
            # coefficient sweeps continuously from today's bot to a learned one.
            model = load_value_model()
            if model is not None:
                total += weights.eval_value_scale * model.advantage(game, owner)
        return total

    @staticmethod
    def _mobile_count(game: Game, side: Owner) -> int:
        """Pieces this side can actually move.

        A side with none of these has lost, whatever material it still holds --
        mines and flags never move, and a headquarters piece is frozen there.
        """
        return sum(
            1
            for position, piece in game.board.items()
            if piece.owner == side
            and piece.kind.movable
            and position not in HEADQUARTERS
        )

    @staticmethod
    def _commander_shield(game: Game, side: Owner) -> float:
        """1 while this side's commander lives and its flag is still hidden.

        The commander's life is worth exactly the concealment it buys, so the
        term vanishes once the flag is out -- at that point the commander is
        just another piece. Evaluated on determinized worlds, so reading the
        opponent's commander is a guess, not a peek.
        """
        squares = game.flag_candidates(side)
        if not squares or any(game.board[square].revealed for square in squares):
            return 0.0
        return float(
            any(
                piece.owner == side and piece.kind == PieceKind.COMMANDER
                for piece in game.board.values()
            )
        )

    def _flag_pressure(
        self, game: Game, owner: Owner, own_moves: list[Move]
    ) -> float:
        """Reward closing on the enemy headquarters, punish losing our own.

        Without this the evaluation is pure material and the bot has no reason
        to ever go for the win.
        """
        weights = self.weights
        value = 0.0
        attacking = game.flag_candidates(owner.other)
        reach = self._closest_raider(game, owner, attacking)
        if reach is not None:
            value += max(0.0, weights.eval_horizon - reach) * (
                weights.eval_hq_attack_certain
                if len(attacking) == 1
                else weights.eval_hq_attack
            )
            if reach <= 1:
                # The mirror of `eval_hq_breach` below: our raider takes the
                # flag next ply unless it is removed, and the linear term above
                # pays the same for closing 6->5 as for closing 2->1. Ships at 0.
                value += weights.eval_hq_storm
        defending = game.flag_candidates(owner)
        threat = self._closest_raider(game, owner.other, defending)
        if threat is not None:
            value -= max(0.0, weights.eval_horizon - threat) * (
                weights.eval_hq_defense_certain
                if len(defending) == 1
                else weights.eval_hq_defense
            )
            if threat <= 1:
                # An enemy standing next to a headquarters that may hold our
                # flag takes it next ply unless we remove it. The linear term
                # above barely distinguishes that from a distant piece.
                value -= weights.eval_hq_breach
            # Measured dead on 2026-08-01: +0.0008 +/- 0.0097 over 806 paired
            # games, an SE tight enough to exclude anything above +0.016. Kept
            # at 0 so the harness can still A/B it (D022), but guarded, because
            # unlike the other retired coefficients this one is not a free
            # multiply -- `_cover` walks the board and the move list on every
            # single evaluation, and the browser's `deep` tier already peaks at
            # 446ms against a 420ms budget.
            if weights.eval_hq_supply and threat <= REINFORCE_RANGE:
                # Supply, not occupancy: a raider this close will be on an
                # approach square within a ply or two, and by then the only
                # thing that matters is how many pieces can still answer.
                cover = self._cover(game, owner, defending, own_moves)
                value -= weights.eval_hq_supply * max(0, MIN_COVER - cover)
        value += weights.eval_hq_guard * self._guards(game, owner, defending)
        return value

    @staticmethod
    def _approaches(squares: list[Position]) -> set[Position]:
        """Squares an attacker must stand on to reach one of ``squares``."""
        return {
            (square[0] + row, square[1] + column)
            for square in squares
            for row, column in ((1, 0), (-1, 0), (0, 1), (0, -1))
            if 0 <= square[0] + row < ROWS and 0 <= square[1] + column < COLUMNS
        }

    @classmethod
    def _cover(
        cls,
        game: Game,
        owner: Owner,
        squares: list[Position],
        own_moves: list[Move],
    ) -> int:
        """Distinct own pieces that hold, or can reach, a flag approach square.

        Counts a piece already standing on an approach *and* a piece one move
        from stepping onto it, because both can answer a raider -- and the
        losses came from having neither. The three-mine screen scores zero
        here: a mine cannot move, so it can never re-take the square it dies
        on.
        """
        if not squares:
            return 0
        approaches = cls._approaches(squares)
        holding = {
            position
            for position, piece in game.board.items()
            if piece.owner == owner
            and piece.kind.movable
            and position not in HEADQUARTERS
            and position in approaches
        }
        arriving = {move.src for move in own_moves if move.dst in approaches}
        return len(holding | arriving)

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

    def _closest_raider(
        self, game: Game, side: Owner, targets: list[Position]
    ) -> int | None:
        """Distance in *moves*, so a raider sitting on a railway counts as near."""
        if not targets:
            return None
        metric = (
            move_distance if self.weights.use_move_distance >= 0.5 else _distance
        )
        return min(
            (
                metric(position, target)
                for position, piece in game.board.items()
                if piece.owner == side
                and piece.kind.movable
                and position not in HEADQUARTERS
                for target in targets
            ),
            default=None,
        )
