from __future__ import annotations

import random
from dataclasses import dataclass, field

from .board import (
    CAMPS,
    HEADQUARTERS,
    engineer_only_move,
    engineer_rail_destinations,
    road_neighbors,
    straight_rail_destinations,
)
from .deployment import headquarters, random_deployment
from .types import Move, Owner, Piece, PieceKind, Position, format_position


@dataclass(frozen=True, slots=True)
class MoveRecord:
    move: Move
    attacker: Piece
    defender: Piece | None
    outcome: int | None
    #: True when no rank but an engineer could have made this move. Stamped in
    #: `apply` because it depends on the board *before* the move, which no later
    #: reader can reconstruct. Public information -- either side may use it.
    engineer_only: bool = False


@dataclass(frozen=True, slots=True)
class ObservedMove:
    move: Move
    attacker_owner: Owner
    had_battle: bool
    outcome: int | None
    own_kind: PieceKind | None
    #: Only an engineer could have made this move. A certainty rather than a
    #: prior, and the one piece of rank information a *quiet* move carries.
    engineer_only: bool = False


@dataclass
class Game:
    board: dict[Position, Piece]
    turn: Owner = Owner.HUMAN
    winner: Owner | None = None
    draw: bool = False
    move_count: int = 0
    history: list[str] = field(default_factory=list)
    records: list[MoveRecord] = field(default_factory=list)

    @classmethod
    def new(
        cls,
        seed: int | None = None,
        first: Owner = Owner.HUMAN,
        bot_deployment: dict[Position, Piece] | None = None,
    ) -> "Game":
        rng = random.Random(seed)
        board: dict[Position, Piece] = {}
        if bot_deployment is None:
            cls._deploy_side(board, Owner.BOT, rng)
        else:
            board.update(bot_deployment)
        cls._deploy_side(board, Owner.HUMAN, rng)
        return cls(board=board, turn=first)

    @staticmethod
    def _deploy_side(
        board: dict[Position, Piece], owner: Owner, rng: random.Random
    ) -> None:
        board.update(random_deployment(owner, rng))

    def clone(self) -> "Game":
        return Game(
            board=self.board.copy(),
            turn=self.turn,
            winner=self.winner,
            draw=self.draw,
            move_count=self.move_count,
            history=self.history.copy(),
            records=self.records.copy(),
        )

    @property
    def over(self) -> bool:
        return self.winner is not None or self.draw

    def legal_moves(self, owner: Owner | None = None) -> list[Move]:
        player = self.turn if owner is None else owner
        if self.over:
            return []
        moves: list[Move] = []
        occupied = set(self.board)
        for source, piece in self.board.items():
            if piece.owner != player or not piece.kind.movable:
                continue
            moves.extend(
                Move(source, target)
                for target in self._destinations(source, occupied)
            )
        return moves

    def _destinations(
        self, source: Position, occupied: set[Position] | None = None
    ) -> set[Position]:
        piece = self.board[source]
        if source in HEADQUARTERS:
            return set()
        if occupied is None:
            occupied = set(self.board)
        destinations = set(road_neighbors(source))
        if piece.kind == PieceKind.ENGINEER:
            destinations |= engineer_rail_destinations(source, occupied)
        else:
            destinations |= straight_rail_destinations(source, occupied)

        valid: set[Position] = set()
        for target in destinations:
            occupant = self.board.get(target)
            if occupant is None:
                valid.add(target)
            elif occupant.owner != piece.owner and target not in CAMPS:
                valid.add(target)
        return valid

    def apply(self, move: Move) -> str:
        if self.over:
            raise ValueError("对局已经结束")
        if move not in self.legal_moves():
            raise ValueError(f"非法走法：{move}")
        # Before the board changes: whether only an engineer could have made
        # this move depends on what was blocking the rails at the time.
        engineer_only = engineer_only_move(move.src, move.dst, set(self.board))
        attacker = self.board.pop(move.src)
        defender = self.board.get(move.dst)
        if defender is None:
            self.board[move.dst] = attacker
            result = f"{move}"
            outcome = None
        else:
            outcome = battle_outcome(attacker.kind, defender.kind)
            result = self._resolve_battle(move, attacker, defender, outcome)

        self.move_count += 1
        self.history.append(result)
        self.records.append(
            MoveRecord(move, attacker, defender, outcome, engineer_only)
        )
        if self.winner is None:
            next_player = self.turn.other
            if not self.legal_moves(next_player):
                self.winner = self.turn
            else:
                self.turn = next_player
        return result

    def _resolve_battle(
        self, move: Move, attacker: Piece, defender: Piece, outcome: int
    ) -> str:
        self.board.pop(move.dst)
        hidden_attacker = Piece(attacker.owner, attacker.kind, False)
        hidden_defender = Piece(defender.owner, defender.kind, False)
        revealed_flags: list[Position] = []
        if attacker.kind == PieceKind.COMMANDER and outcome <= 0:
            flag = self._reveal_flag(attacker.owner)
            if flag is not None:
                revealed_flags.append(flag)
        if defender.kind == PieceKind.COMMANDER and outcome >= 0:
            flag = self._reveal_flag(defender.owner)
            if flag is not None:
                revealed_flags.append(flag)
        if defender.kind == PieceKind.FLAG:
            self.winner = attacker.owner
        elif outcome > 0:
            self.board[move.dst] = hidden_attacker
        elif outcome < 0:
            self.board[move.dst] = hidden_defender
        if defender.kind == PieceKind.FLAG:
            return f"{move}: 军旗被夺，对局结束"
        result = "攻方获胜" if outcome > 0 else "守方获胜" if outcome < 0 else "同归于尽"
        message = f"{move}: {result}（双方军衔保持隐藏）"
        if revealed_flags:
            locations = "、".join(format_position(position) for position in revealed_flags)
            message += f"；司令阵亡，军旗位置亮出：{locations}"
        return message

    def observations(self, owner: Owner, after: int = 0) -> list[ObservedMove]:
        observations: list[ObservedMove] = []
        for record in self.records[after:]:
            own_kind: PieceKind | None = None
            if record.attacker.owner == owner:
                own_kind = record.attacker.kind
            elif record.defender is not None and record.defender.owner == owner:
                own_kind = record.defender.kind
            observations.append(
                ObservedMove(
                    move=record.move,
                    attacker_owner=record.attacker.owner,
                    had_battle=record.defender is not None,
                    outcome=record.outcome,
                    own_kind=own_kind,
                    engineer_only=record.engineer_only,
                )
            )
        return observations

    def _reveal_flag(self, owner: Owner) -> Position | None:
        position = self.flag_position(owner)
        if position is None:
            return None
        flag = self.board[position]
        self.board[position] = Piece(flag.owner, flag.kind, True)
        return position

    def flag_position(self, owner: Owner) -> Position | None:
        for position, piece in self.board.items():
            if piece.owner == owner and piece.kind == PieceKind.FLAG:
                return position
        return None

    def flag_candidates(self, owner: Owner) -> list[Position]:
        """Headquarters squares that may still hold ``owner``'s flag.

        Deployment fills every non-camp home square, headquarters pieces can
        never move, and a flag never leaves its headquarters -- so while the
        game is running the flag is under one of the own-side pieces still
        standing on a headquarters square. Reads occupancy only, never a rank,
        so either side may use it against the other.
        """
        return sorted(
            position
            for position in headquarters(owner)
            if (piece := self.board.get(position)) is not None and piece.owner == owner
        )


def battle_outcome(attacker: PieceKind, defender: PieceKind) -> int:
    """Return positive for attacker win, negative for defender win, zero for both."""
    if defender == PieceKind.FLAG:
        return 1
    if attacker == PieceKind.BOMB or defender == PieceKind.BOMB:
        return 0
    if defender == PieceKind.MINE:
        return 1 if attacker == PieceKind.ENGINEER else -1
    if attacker == defender:
        return 0
    return 1 if attacker.value < defender.value else -1
