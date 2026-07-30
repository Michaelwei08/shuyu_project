from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum


class Owner(IntEnum):
    HUMAN = 0
    BOT = 1

    @property
    def other(self) -> "Owner":
        return Owner.BOT if self == Owner.HUMAN else Owner.HUMAN


class PieceKind(IntEnum):
    FLAG = 0
    COMMANDER = 1
    GENERAL = 2
    MAJOR_GENERAL = 3
    BRIGADIER = 4
    COLONEL = 5
    MAJOR = 6
    CAPTAIN = 7
    LIEUTENANT = 8
    ENGINEER = 9
    MINE = 10
    BOMB = 11

    @property
    def movable(self) -> bool:
        return self not in (PieceKind.FLAG, PieceKind.MINE)


PIECE_COUNTS: dict[PieceKind, int] = {
    PieceKind.FLAG: 1,
    PieceKind.COMMANDER: 1,
    PieceKind.GENERAL: 1,
    PieceKind.MAJOR_GENERAL: 2,
    PieceKind.BRIGADIER: 2,
    PieceKind.COLONEL: 2,
    PieceKind.MAJOR: 2,
    PieceKind.CAPTAIN: 3,
    PieceKind.LIEUTENANT: 3,
    PieceKind.ENGINEER: 3,
    PieceKind.MINE: 3,
    PieceKind.BOMB: 2,
}

SYMBOLS: dict[PieceKind, str] = {
    PieceKind.FLAG: "旗",
    PieceKind.COMMANDER: "司",
    PieceKind.GENERAL: "军",
    PieceKind.MAJOR_GENERAL: "师",
    PieceKind.BRIGADIER: "旅",
    PieceKind.COLONEL: "团",
    PieceKind.MAJOR: "营",
    PieceKind.CAPTAIN: "连",
    PieceKind.LIEUTENANT: "排",
    PieceKind.ENGINEER: "工",
    PieceKind.MINE: "雷",
    PieceKind.BOMB: "炸",
}


@dataclass(frozen=True, slots=True)
class Piece:
    owner: Owner
    kind: PieceKind
    revealed: bool = False


Position = tuple[int, int]


@dataclass(frozen=True, slots=True)
class Move:
    src: Position
    dst: Position

    def __str__(self) -> str:
        return f"{format_position(self.src)}-{format_position(self.dst)}"


def format_position(position: Position) -> str:
    row, column = position
    return f"{chr(ord('A') + column)}{row + 1}"


def parse_position(text: str) -> Position:
    value = text.strip().upper()
    if len(value) < 2 or value[0] not in "ABCDE" or not value[1:].isdigit():
        raise ValueError(f"无效坐标：{text!r}，请使用 A1 到 E12")
    position = (int(value[1:]) - 1, ord(value[0]) - ord("A"))
    if not (0 <= position[0] < 12):
        raise ValueError(f"无效坐标：{text!r}，行号必须为 1 到 12")
    return position


def parse_move(text: str) -> Move:
    normalized = text.replace("→", "-").replace(" ", "-")
    parts = [part for part in normalized.split("-") if part]
    if len(parts) != 2:
        raise ValueError("走法格式应为 A10-A9")
    return Move(parse_position(parts[0]), parse_position(parts[1]))
