from __future__ import annotations

import json
import random
from collections import Counter
from pathlib import Path

from .board import CAMPS
from .types import (
    Owner,
    PIECE_COUNTS,
    Piece,
    PieceKind,
    Position,
    format_position,
    parse_position,
)


def deployment_positions(owner: Owner) -> list[Position]:
    rows = range(0, 6) if owner == Owner.BOT else range(6, 12)
    return [
        (row, column)
        for row in rows
        for column in range(5)
        if (row, column) not in CAMPS
    ]


def headquarters(owner: Owner) -> set[Position]:
    return {(0, 1), (0, 3)} if owner == Owner.BOT else {(11, 1), (11, 3)}


def rear_rows(owner: Owner) -> set[int]:
    return {0, 1} if owner == Owner.BOT else {10, 11}


def front_row(owner: Owner) -> int:
    return 5 if owner == Owner.BOT else 6


def validate_deployment(
    pieces: dict[Position, Piece], owner: Owner
) -> list[str]:
    errors: list[str] = []
    allowed = set(deployment_positions(owner))
    own_pieces = {
        position: piece for position, piece in pieces.items() if piece.owner == owner
    }
    if set(own_pieces) - allowed:
        errors.append("棋子只能放在己方兵站，不能放在行营或对方区域")

    counts = Counter(piece.kind for piece in own_pieces.values())
    if counts != Counter(PIECE_COUNTS):
        errors.append("棋子数量或编制不正确")

    flag_positions = [
        position
        for position, piece in own_pieces.items()
        if piece.kind == PieceKind.FLAG
    ]
    if any(position not in headquarters(owner) for position in flag_positions):
        errors.append("军旗必须放在己方两个大本营之一")

    mine_positions = [
        position
        for position, piece in own_pieces.items()
        if piece.kind == PieceKind.MINE
    ]
    if any(position[0] not in rear_rows(owner) for position in mine_positions):
        errors.append("地雷只能放在最后两排")

    bomb_positions = [
        position
        for position, piece in own_pieces.items()
        if piece.kind == PieceKind.BOMB
    ]
    if any(position[0] == front_row(owner) for position in bomb_positions):
        errors.append("炸弹不能放在最前排")
    return errors


def random_deployment(owner: Owner, rng: random.Random) -> dict[Position, Piece]:
    available = deployment_positions(owner)
    result: dict[Position, Piece] = {}

    flag_position = rng.choice(sorted(headquarters(owner)))
    result[flag_position] = Piece(owner, PieceKind.FLAG)
    available.remove(flag_position)

    mine_slots = [position for position in available if position[0] in rear_rows(owner)]
    for position in rng.sample(mine_slots, PIECE_COUNTS[PieceKind.MINE]):
        result[position] = Piece(owner, PieceKind.MINE)
        available.remove(position)

    bomb_slots = [position for position in available if position[0] != front_row(owner)]
    for position in rng.sample(bomb_slots, PIECE_COUNTS[PieceKind.BOMB]):
        result[position] = Piece(owner, PieceKind.BOMB)
        available.remove(position)

    remaining: list[PieceKind] = []
    for kind, count in PIECE_COUNTS.items():
        if kind not in (PieceKind.FLAG, PieceKind.MINE, PieceKind.BOMB):
            remaining.extend([kind] * count)
    rng.shuffle(remaining)
    rng.shuffle(available)
    for position, kind in zip(available, remaining, strict=True):
        result[position] = Piece(owner, kind)

    errors = validate_deployment(result, owner)
    if errors:
        raise RuntimeError("随机布阵生成失败：" + "；".join(errors))
    return result


#: Most mines allowed on the flag headquarters' own neighbours. Was effectively
#: 3 (a full seal) until ten replayed games showed the seal losing every time;
#: see `_build_strategic`. Kept as a parameter so the paired harness can A/B a
#: capped screen against the old full seal.
SCREEN_MINE_CAP = 2


def strategic_deployment(
    owner: Owner,
    rng: random.Random,
    attempts: int = 60,
    screen_cap: int | None = None,
) -> dict[Position, Piece]:
    """A fresh, legal, non-random-looking opening.

    A fixed opening is worth nothing once the opponent has seen it twice, and a
    uniform random one wastes material. This keeps the shape sensible -- mines
    screening the flag, a cheap decoy in the unused headquarters, leaders off
    the back rows -- while varying every game.

    ``screen_cap`` overrides :data:`SCREEN_MINE_CAP`; pass 3 to reproduce the
    old full seal.
    """
    for _ in range(attempts):
        result = _build_strategic(owner, rng, screen_cap)
        if result is not None and not validate_deployment(result, owner):
            return result
    return random_deployment(owner, rng)


def _build_strategic(
    owner: Owner, rng: random.Random, screen_cap: int | None = None
) -> dict[Position, Piece] | None:
    result: dict[Position, Piece] = {}
    free = set(deployment_positions(owner))
    rear = rear_rows(owner)

    def place(position: Position, kind: PieceKind) -> None:
        result[position] = Piece(owner, kind)
        free.discard(position)

    flag_hq, decoy_hq = rng.sample(sorted(headquarters(owner)), 2)
    place(flag_hq, PieceKind.FLAG)
    # Whatever sits in the other headquarters is frozen there for the whole
    # game, so it must be a piece we can afford to lose.
    place(decoy_hq, rng.choice([PieceKind.LIEUTENANT, PieceKind.CAPTAIN]))

    rear_free = [position for position in free if position[0] in rear]
    screen = [
        position
        for position in rear_free
        if abs(position[0] - flag_hq[0]) + abs(position[1] - flag_hq[1]) == 1
    ]
    rng.shuffle(screen)
    rng.shuffle(rear_free)
    # A headquarters has exactly three orthogonal neighbours, all of them in the
    # rear rows, so all three mines *can* seal the flag. Ten replayed games say
    # not to. Those neighbours are alternative doors, not three locks on one
    # door -- the attacker only has to open the cheapest, and in eight of eight
    # flag losses exactly one mine was cleared, always the square the killer
    # then stood on. Sealing also costs more than it buys: a mine cannot move,
    # so a full screen leaves no square from which the bot can ever post a
    # mobile defender, and `eval_hq_guard` is pinned at zero until the bot's own
    # screen has been destroyed. Leave at least one neighbour to a piece that
    # can fight back and be replaced.
    cap = SCREEN_MINE_CAP if screen_cap is None else screen_cap
    guards = screen[: cap if rng.random() < 0.75 else max(1, cap - 1)]
    # Exclude the whole screen from the tail, not just the chosen guards, or the
    # remaining mines land back on the neighbours this cap exists to keep free.
    mine_slots = guards + [
        position for position in rear_free if position not in screen
    ]
    if len(mine_slots) < PIECE_COUNTS[PieceKind.MINE]:
        return None
    for position in mine_slots[: PIECE_COUNTS[PieceKind.MINE]]:
        place(position, PieceKind.MINE)

    # Bombs are wasted on the back rank and illegal on the front one.
    midfield = [
        position
        for position in free
        if position[0] != front_row(owner) and position[0] not in rear
    ]
    rng.shuffle(midfield)
    if len(midfield) < PIECE_COUNTS[PieceKind.BOMB] + 2:
        return None
    for position in midfield[: PIECE_COUNTS[PieceKind.BOMB]]:
        place(position, PieceKind.BOMB)

    # Keep the commander and general out of the rear rows so they can be used.
    leaders = [PieceKind.COMMANDER, PieceKind.GENERAL]
    forward = [position for position in free if position[0] not in rear]
    rng.shuffle(forward)
    if len(forward) < len(leaders):
        return None
    for position, kind in zip(forward, leaders, strict=False):
        place(position, kind)

    placed = Counter(piece.kind for piece in result.values())
    remaining: list[PieceKind] = []
    for kind, count in PIECE_COUNTS.items():
        remaining.extend([kind] * (count - placed[kind]))
    leftover = sorted(free)
    if len(leftover) != len(remaining):
        return None
    rng.shuffle(remaining)
    for position, kind in zip(leftover, remaining, strict=True):
        place(position, kind)
    return result


def swap_pieces(
    board: dict[Position, Piece],
    owner: Owner,
    left: Position,
    right: Position,
) -> None:
    if left == right:
        raise ValueError("请选择两个不同的位置")
    if left not in board or right not in board:
        raise ValueError("两个位置都必须有棋子")
    if board[left].owner != owner or board[right].owner != owner:
        raise ValueError("只能交换自己的棋子")
    board[left], board[right] = board[right], board[left]
    errors = validate_deployment(board, owner)
    if errors:
        board[left], board[right] = board[right], board[left]
        raise ValueError("该交换会产生非法阵型：" + "；".join(errors))


def save_deployment(pieces: dict[Position, Piece], path: str | Path) -> None:
    errors = validate_deployment(pieces, Owner.BOT)
    if errors:
        raise ValueError("不能保存非法 Bot 阵型：" + "；".join(errors))
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        format_position(position): piece.kind.name
        for position, piece in sorted(pieces.items())
        if piece.owner == Owner.BOT
    }
    with destination.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)


def load_deployment(path: str | Path) -> dict[Position, Piece]:
    with Path(path).open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    pieces = {
        parse_position(position): Piece(Owner.BOT, PieceKind[kind])
        for position, kind in payload.items()
    }
    errors = validate_deployment(pieces, Owner.BOT)
    if errors:
        raise ValueError("Bot 布局模型不合法：" + "；".join(errors))
    return pieces
