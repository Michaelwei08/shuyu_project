from __future__ import annotations

import json
import random
from collections import Counter
from dataclasses import dataclass
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


#: Most mines allowed on the flag headquarters' own neighbours. Three is a full
#: seal. Dropping this to 2 was tried on 2026-08-01 and **measured much worse**
#: -- see `_build_strategic`. Kept as a parameter so the paired harness can A/B
#: it; that is how the attempt was caught.
SCREEN_MINE_CAP = 3


@dataclass(frozen=True)
class DeploymentFamily:
    """One point in the opening's strategy space.

    Deployment is a *simultaneous one-shot* game played before a move is made,
    so its solution is a distribution over openings rather than a single best
    one. Everything the generator used to hardcode lives here instead, so a
    family can be named on the command line, carried in a picklable `Job`, and
    mixed over -- see :mod:`junqi.deployment_game`.
    """

    name: str
    #: Kinds allowed in the headquarters that does *not* hold the flag. Whatever
    #: goes there is frozen for the whole game and its square is public, so this
    #: choice sets the price of probing a headquarters.
    decoy: tuple[PieceKind, ...] = (PieceKind.LIEUTENANT, PieceKind.CAPTAIN)
    screen_cap: int = SCREEN_MINE_CAP
    #: Bombs placed on the flag headquarters' own neighbours, taken out of the
    #: midfield allocation. A bomb on a door kills any attacker once; a mine
    #: kills every attacker except an engineer, forever.
    home_bombs: int = 0
    #: Ignore all of the above and deploy uniformly at random. The floor, and
    #: the control that says how much the shaped generator is worth at all.
    uniform: bool = False
    #: A distribution over other families, drawn once per game. This is the
    #: point of the exercise: the opening is chosen simultaneously and never
    #: seen, so its solution is a mixed strategy, and a pure family is only the
    #: degenerate case. Drawn from the deploying side's own stream, so mixing
    #: one army never perturbs the other.
    mixture: tuple[tuple[str, float], ...] = ()
    #: Mine the full screen every game, instead of dropping to `cap - 1` a
    #: quarter of the time. **Now the default, on measurement**: the generator
    #: used to seal all three doors only 75% of the time, leaving a mobile piece
    #: on a flag door in 26.7% of openings -- and a mobile guard is one the
    #: search walks off for any capture worth more than `eval_hq_guard = 5.5`.
    #: It was caught doing exactly that in a human game on 2026-08-04: a major
    #: general left the door at ply 14 to take a 3-point engineer, and the flag
    #: fell through that square 27 plies later.
    #:
    #: Sealing unconditionally measured **+0.0384 +/- 0.0067 (clustered 0.0089),
    #: p ~ 0, over 2400 paired games** -- 72.8% against 69.6%, with the held-out
    #: four echoing it at +0.0365. `seal-75` is the old behaviour, kept as the
    #: switched-off end.
    always_seal: bool = True


FAMILIES: dict[str, DeploymentFamily] = {
    # The shipped generator, as of 2026-08-02.
    "standard": DeploymentFamily("standard"),
    # A probe of the decoy headquarters currently kills a lieutenant or captain
    # and buys certainty about where the flag is. These make it cost something.
    "decoy-mine": DeploymentFamily("decoy-mine", decoy=(PieceKind.MINE,)),
    "decoy-bomb": DeploymentFamily("decoy-bomb", decoy=(PieceKind.BOMB,)),
    # Known bad: measured -0.0985 +/- 0.0217 over 806 paired games. Kept in the
    # family set on purpose -- a matrix that cannot rediscover a known-bad
    # strategy is not measuring anything.
    "screen2": DeploymentFamily("screen2", screen_cap=2),
    "bomb-home": DeploymentFamily("bomb-home", home_bombs=1),
    # The generator as it stood before 2026-08-04: all three doors sealed only
    # 75% of the time, so 26.7% of openings hand a flag door to a mobile piece.
    # Kept as the off end of the switch that replaced it, the same way
    # `unknown_risk` and `screen2` are kept -- a superseded variant that cannot
    # be selected is a result nobody can re-check.
    "seal-75": DeploymentFamily("seal-75", always_seal=False),
    # Identical to `standard` now that sealing is the default. Retained because
    # the command that measured the change still names it, and because an
    # identical family is a free instrument check: it must score exactly
    # +0.0000 +/- 0.0000 against `standard`.
    "seal-always": DeploymentFamily("seal-always", always_seal=True),
    "uniform": DeploymentFamily("uniform", uniform=True),
}

DEFAULT_FAMILY = "standard"


def parse_mixture(text: str) -> DeploymentFamily:
    """``mix:standard=0.4,decoy-mine=0.6`` -> a family that draws per game.

    Spelled inline rather than loaded from a model file so that the exact
    distribution under test is visible in the command line, in the log and in
    the pickled `Job` -- a mixture read from disk inside a worker is a silent
    dependency on whatever that file happened to contain at the time.
    """
    parts = []
    for item in text[len("mix:") :].split(","):
        name, _, share = item.partition("=")
        if name not in FAMILIES:
            raise KeyError(f"未知布阵族：{name}（可选 {sorted(FAMILIES)}）")
        parts.append((name, float(share)))
    total = sum(share for _, share in parts)
    if total <= 0:
        raise ValueError(f"混合权重之和必须为正：{text}")
    return DeploymentFamily(
        text, mixture=tuple((name, share / total) for name, share in parts)
    )


def resolve_family(family: str | int | None) -> DeploymentFamily:
    """Accept a family name, an inline mixture, a bare screen cap, or nothing.

    The bare int is the older `--screen-cap` interface, which predates families
    and is still what `models/ab/` comparisons and the existing tests speak.
    """
    if family is None:
        return FAMILIES[DEFAULT_FAMILY]
    if isinstance(family, int):
        # `always_seal=False` on purpose: `--screen-cap` predates the 2026-08-04
        # change, and the -0.0985 cap-2-versus-cap-3 result was measured with the
        # 75% branch in place. Keeping the old semantics here keeps that number
        # reproducible instead of silently re-defining what it compared.
        return DeploymentFamily(
            f"screen{family}", screen_cap=family, always_seal=False
        )
    if family.startswith("mix:"):
        return parse_mixture(family)
    if family not in FAMILIES:
        raise KeyError(f"未知布阵族：{family}（可选 {sorted(FAMILIES)}）")
    return FAMILIES[family]


def strategic_deployment(
    owner: Owner,
    rng: random.Random,
    attempts: int = 60,
    screen_cap: int | None = None,
    family: str | int | None = None,
) -> dict[Position, Piece]:
    """A fresh, legal, non-random-looking opening.

    A fixed opening is worth nothing once the opponent has seen it twice, and a
    uniform random one wastes material. This keeps the shape sensible -- mines
    screening the flag, a cheap decoy in the unused headquarters, leaders off
    the back rows -- while varying every game.

    ``screen_cap`` overrides :data:`SCREEN_MINE_CAP`; pass 3 to reproduce the
    old full seal. ``family`` names a whole :class:`DeploymentFamily` and takes
    precedence over it.
    """
    chosen = resolve_family(family if family is not None else screen_cap)
    if chosen.mixture:
        draw = rng.random()
        for name, share in chosen.mixture:
            draw -= share
            if draw <= 0:
                chosen = FAMILIES[name]
                break
        else:  # pragma: no cover - float dust on the last component
            chosen = FAMILIES[chosen.mixture[-1][0]]
    if chosen.uniform:
        return random_deployment(owner, rng)
    for _ in range(attempts):
        result = _build_strategic(owner, rng, chosen)
        if result is not None and not validate_deployment(result, owner):
            return result
    return random_deployment(owner, rng)


def _build_strategic(
    owner: Owner, rng: random.Random, family: DeploymentFamily | None = None
) -> dict[Position, Piece] | None:
    chosen = family if family is not None else FAMILIES[DEFAULT_FAMILY]
    result: dict[Position, Piece] = {}
    free = set(deployment_positions(owner))
    rear = rear_rows(owner)

    def place(position: Position, kind: PieceKind) -> None:
        result[position] = Piece(owner, kind)
        free.discard(position)

    def budget(kind: PieceKind) -> int:
        placed = sum(1 for piece in result.values() if piece.kind == kind)
        return PIECE_COUNTS[kind] - placed

    flag_hq, decoy_hq = rng.sample(sorted(headquarters(owner)), 2)
    place(flag_hq, PieceKind.FLAG)
    # Whatever sits in the other headquarters is frozen there for the whole
    # game, so it must be a piece we can afford to lose -- or one that makes
    # the probe expensive, which is what the non-default families try.
    place(decoy_hq, rng.choice(list(chosen.decoy)))

    rear_free = [position for position in free if position[0] in rear]
    screen = [
        position
        for position in rear_free
        if abs(position[0] - flag_hq[0]) + abs(position[1] - flag_hq[1]) == 1
    ]
    rng.shuffle(screen)
    rng.shuffle(rear_free)
    # A headquarters has exactly three orthogonal neighbours, all of them in the
    # rear rows, so all three mines can seal the flag.
    #
    # Capping this at 2 looked obviously right and measured **-0.0985 +/-
    # 0.0217 over 806 paired games** -- about eight points of win rate, the
    # largest single effect in that round and four times the size of anything
    # it was bundled with. The argument for capping was that the neighbours are
    # alternative doors rather than three locks on one door, since in eight of
    # eight replayed flag losses exactly one mine was cleared, always the square
    # the killer then stood on. That observation is true and the inference from
    # it was wrong: instrumenting 208 pool games showed the flag falling at
    # essentially the same rate under both caps (23/104 capped vs 21/104
    # sealed), so the loss was never about flag defence. A mine on a door also
    # never leaves; an ordinary piece on a door is one the bot walks away from
    # for any capture worth more than `eval_hq_guard`.
    cap = chosen.screen_cap
    # A bomb on a door takes one of the doors the mines would otherwise seal,
    # and comes out of the midfield allocation rather than adding material.
    for position in screen[: min(chosen.home_bombs, budget(PieceKind.BOMB))]:
        place(position, PieceKind.BOMB)
    screen = [position for position in screen if position in free]
    # The 0.75 is what leaves a quarter of openings with a walkable door. Drawn
    # unconditionally so `always_seal` consumes the same number of values from
    # the stream as `standard` does -- otherwise the two families would produce
    # different *layouts* from one seed for a second reason, and the comparison
    # would not be measuring the seal.
    loosen = rng.random() >= 0.75
    guards = screen[: max(1, cap - 1) if loosen and not chosen.always_seal else cap]
    # Exclude the whole screen from the tail, not just the chosen guards, or the
    # remaining mines land back on the neighbours this cap exists to keep free.
    mine_slots = guards + [
        position for position in rear_free if position not in screen and position in free
    ]
    # A family that spends a mine on the decoy headquarters has fewer left to
    # place here, which is exactly the trade it is making.
    mines = budget(PieceKind.MINE)
    if len(mine_slots) < mines:
        return None
    for position in mine_slots[:mines]:
        place(position, PieceKind.MINE)

    # Bombs are wasted on the back rank and illegal on the front one.
    midfield = [
        position
        for position in free
        if position[0] != front_row(owner) and position[0] not in rear
    ]
    rng.shuffle(midfield)
    bombs = budget(PieceKind.BOMB)
    if len(midfield) < bombs + 2:
        return None
    for position in midfield[:bombs]:
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


def save_deployment(
    pieces: dict[Position, Piece], path: str | Path, owner: Owner = Owner.BOT
) -> None:
    """Write one side's layout to disk.

    ``owner`` defaults to the bot so `layout_training` and `--deployment-model`
    keep their historical meaning, but a person's own opening is the same object
    on the other half of the board -- worth saving once rather than re-swapping
    it by hand every game.
    """
    errors = validate_deployment(pieces, owner)
    if errors:
        raise ValueError(f"不能保存非法阵型（{owner.name}）：" + "；".join(errors))
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        format_position(position): piece.kind.name
        for position, piece in sorted(pieces.items())
        if piece.owner == owner
    }
    with destination.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)


def load_deployment(
    path: str | Path, owner: Owner = Owner.BOT
) -> dict[Position, Piece]:
    """Read a layout saved by :func:`save_deployment`.

    The squares in the file already imply which half of the board it is for, so
    a bot layout loaded as a human one fails `validate_deployment` rather than
    silently producing a board with 50 pieces on one side.
    """
    with Path(path).open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    pieces = {
        parse_position(position): Piece(owner, PieceKind[kind])
        for position, kind in payload.items()
    }
    errors = validate_deployment(pieces, owner)
    if errors:
        raise ValueError(
            f"{path} 不是合法的 {owner.name} 阵型：" + "；".join(errors)
        )
    return pieces
