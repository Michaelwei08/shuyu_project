"""The one place an LLM-facing prompt is allowed to read a board.

Everything an agent is told about a position must come from :class:`Observation`,
and :func:`build_observation` is the only function that touches ``game.board``.
That is the whole point of this module: the engine's central invariant is that a
bot never reads a hidden enemy rank, and a text serialiser is by far the easiest
place to break it by accident.

Two landmines, both of which look harmless:

* ``game.legal_moves(enemy)`` is **not** public information. It filters on
  ``piece.kind.movable``, so its output tells you which enemy squares hold a
  mine or the flag. Only ever call it for ``me``.
* ``Piece.revealed`` *is* public -- it is set only by the commander-death flag
  reveal -- and a revealed piece is therefore always a flag. That is the single
  enemy rank a prompt may name.

:class:`EnemySquare` carries no ``kind`` field at all, which is what makes the
guarantee structural rather than a matter of remembering. The paired test is
``test_the_prompt_is_invariant_under_permuting_hidden_ranks``: relabel every
unrevealed enemy piece and the rendered text must come back byte-identical.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from junqi.board import CAMPS, HEADQUARTERS, RAILWAYS
from junqi.game import Game
from junqi.knowledge import MOVABLE_KINDS
from junqi.types import (
    SYMBOLS,
    Move,
    Owner,
    PieceKind,
    Position,
    format_position,
    parse_move,
    parse_position,
)

#: Rows each side deploys on, in display coordinates (1-based).
SIDE_ROWS: dict[Owner, str] = {Owner.HUMAN: "7-12", Owner.BOT: "1-6"}
SIDE_NAMES: dict[Owner, str] = {Owner.HUMAN: "南方", Owner.BOT: "北方"}

#: How many ranks an enemy piece could be before anything has been deduced.
#: Taken from the engine rather than restated, so it tracks `PieceKind.movable`.
MOVABLE_KINDS_COUNT = len(MOVABLE_KINDS)


@dataclass(frozen=True, slots=True)
class OwnPiece:
    """One of your own pieces. You know your own army, so the rank is here."""

    square: Position
    kind: PieceKind


@dataclass(frozen=True, slots=True)
class EnemySquare:
    """An enemy piece as it is legally visible: a location, and nothing else.

    There is deliberately no ``kind`` field. ``revealed`` is safe to carry
    because only a flag is ever revealed.
    """

    square: Position
    revealed: bool


@dataclass(frozen=True, slots=True)
class LoggedMove:
    """One entry of the public move log, from ``Game.observations``.

    ``outcome`` is signed from the *attacker's* view (>0 attacker won, <0
    defender won, 0 mutual destruction) and ``own_kind`` is filled in only when
    the piece involved was yours -- exactly the channel the engine's own bots
    read.
    """

    ply: int
    src: Position
    dst: Position
    by_me: bool
    had_battle: bool
    outcome: int | None
    own_kind: PieceKind | None


@dataclass(frozen=True, slots=True)
class Observation:
    """Everything one player may legally know at one ply."""

    me: Owner
    ply: int
    own: tuple[OwnPiece, ...]
    enemy: tuple[EnemySquare, ...]
    log: tuple[LoggedMove, ...]
    legal: tuple[Move, ...]
    own_flag_squares: tuple[Position, ...]
    enemy_flag_squares: tuple[Position, ...]
    #: Ranks each enemy square could still hold, for the squares where anything
    #: has actually been deduced. Empty when nothing has been.
    belief: tuple[tuple[Position, tuple[PieceKind, ...]], ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {
            "me": self.me.name,
            "ply": self.ply,
            "own": [
                {"square": format_position(item.square), "kind": item.kind.name}
                for item in self.own
            ],
            "enemy": [
                {"square": format_position(item.square), "revealed": item.revealed}
                for item in self.enemy
            ],
            "log": [
                {
                    "ply": item.ply,
                    "move": f"{format_position(item.src)}-{format_position(item.dst)}",
                    "by_me": item.by_me,
                    "had_battle": item.had_battle,
                    "outcome": item.outcome,
                    "own_kind": None if item.own_kind is None else item.own_kind.name,
                }
                for item in self.log
            ],
            "legal": [str(move) for move in self.legal],
            "own_flag_squares": [format_position(s) for s in self.own_flag_squares],
            "enemy_flag_squares": [format_position(s) for s in self.enemy_flag_squares],
            "belief": [
                {
                    "square": format_position(square),
                    "kinds": [kind.name for kind in kinds],
                }
                for square, kinds in self.belief
            ],
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "Observation":
        return cls(
            me=Owner[data["me"]],
            ply=int(data["ply"]),
            own=tuple(
                OwnPiece(parse_position(item["square"]), PieceKind[item["kind"]])
                for item in data["own"]
            ),
            enemy=tuple(
                EnemySquare(parse_position(item["square"]), bool(item["revealed"]))
                for item in data["enemy"]
            ),
            log=tuple(
                LoggedMove(
                    ply=int(item["ply"]),
                    src=parse_move(item["move"]).src,
                    dst=parse_move(item["move"]).dst,
                    by_me=bool(item["by_me"]),
                    had_battle=bool(item["had_battle"]),
                    outcome=item["outcome"],
                    own_kind=(
                        None
                        if item["own_kind"] is None
                        else PieceKind[item["own_kind"]]
                    ),
                )
                for item in data["log"]
            ),
            legal=tuple(parse_move(text) for text in data["legal"]),
            own_flag_squares=tuple(
                parse_position(s) for s in data["own_flag_squares"]
            ),
            enemy_flag_squares=tuple(
                parse_position(s) for s in data["enemy_flag_squares"]
            ),
            belief=tuple(
                (
                    parse_position(item["square"]),
                    tuple(PieceKind[name] for name in item["kinds"]),
                )
                for item in data.get("belief", ())
            ),
        )


def build_observation(
    game: Game,
    me: Owner,
    belief: dict[Position, frozenset[PieceKind]] | None = None,
) -> Observation:
    """Project ``game`` down to what ``me`` may legally know.

    The only function in this package that reads ``game.board``. Enemy entries
    keep their square and their ``revealed`` bit; the rank is dropped on the
    floor here rather than filtered out downstream.
    """
    own: list[OwnPiece] = []
    enemy: list[EnemySquare] = []
    for square, piece in sorted(game.board.items()):
        if piece.owner == me:
            own.append(OwnPiece(square, piece.kind))
        else:
            enemy.append(EnemySquare(square, piece.revealed))

    log = [
        LoggedMove(
            ply=index,
            src=event.move.src,
            dst=event.move.dst,
            by_me=event.attacker_owner == me,
            had_battle=event.had_battle,
            outcome=event.outcome,
            own_kind=event.own_kind,
        )
        for index, event in enumerate(game.observations(me))
    ]

    live = {square for square, _ in ((p.square, p) for p in enemy)}
    deduced = tuple(
        (square, tuple(sorted(kinds, key=lambda kind: kind.value)))
        for square, kinds in sorted((belief or {}).items())
        if square in live
    )

    return Observation(
        me=me,
        ply=game.move_count,
        own=tuple(own),
        enemy=tuple(enemy),
        log=tuple(log),
        # Never `legal_moves(me.other)` -- that leaks which squares are immovable.
        legal=tuple(sorted(game.legal_moves(me), key=str)),
        own_flag_squares=tuple(game.flag_candidates(me)),
        enemy_flag_squares=tuple(game.flag_candidates(me.other)),
        belief=deduced,
    )


@dataclass(frozen=True, slots=True)
class Scaffold:
    """Which legally-available facts get spelled out in the prompt.

    Every field is information the agent already holds implicitly; turning one
    on saves it a deduction rather than telling it anything new. That is what
    makes this a clean scaffolding ladder instead of an information ladder.
    """

    name: str
    #: Enumerate every legal move.
    legal_moves: bool = False
    #: State the headquarters squares the enemy flag can still be under.
    flag_candidates: bool = False
    #: State the deduced rank sets for enemy squares.
    belief: bool = False
    #: Restate your own army as a list as well as drawing it on the grid.
    piece_list: bool = True
    #: Include the rules preamble. Kept as a knob so it can be ablated.
    rules: bool = True


SCAFFOLDS: dict[str, Scaffold] = {
    "raw": Scaffold("raw"),
    "legal": Scaffold("legal", legal_moves=True),
    "derived": Scaffold(
        "derived", legal_moves=True, flag_candidates=True, belief=True
    ),
}


# Placed first in every prompt and byte-identical across calls, so a provider's
# prompt cache can cover it. Do not interpolate anything into this string.
RULES = """\
你在下一盘 12 行 5 列的军棋。列为 A-E，行为 1-12，坐标写作 B7。

棋子与胜负：
- 军衔从大到小：司令、军长、师长、旅长、团长、营长、连长、排长、工兵。军衔大的吃小的。
- 同级别相撞，同归于尽。
- 地雷不能移动。除工兵外，任何棋子碰地雷都会被炸死；只有工兵能排雷。
- 炸弹与任何棋子相撞都同归于尽。
- 军旗不能移动。军旗被吃，该方立即输棋。

棋盘固定格局（下面的格子在图上被棋子占住时不会另外标出，请以这份清单为准）：
- 大本营共 4 格：B1、D1（北方），B12、D12（南方）。
- 行营共 10 格：B3、D3、C4、B5、D5（北方），B8、D8、C9、B10、D10（南方）。
- 铁路格：第 2、6、7、11 行的全部 5 格；A 列和 E 列的第 2 到第 11 行；外加 C6、C7。
  其余格子都是公路。

棋盘规则：
- 行营是安全格：可以走进空行营，但绝不能吃行营里的棋子。
- 大本营里的棋子永远不能移动。军旗一定在本方两个大本营之一。
- 公路上，任何棋子每步只能走到上下左右相邻的一格；行营与斜向相邻的格子之间也可以走。
- 铁路上，除工兵外的棋子只能沿铁路直线滑行，不能拐弯，且不能越过任何棋子。
- 工兵在铁路上可以任意拐弯，同样不能越过棋子。
- 河界在第 6 行与第 7 行之间，只有 A、C、E 三列可以通过。

情报规则（最重要）：
- 双方棋子的军衔始终互相隐藏。
- 战斗只公布结果：攻方获胜、守方获胜、或同归于尽。存活者的军衔仍然不公开。
- 唯一例外：某方司令阵亡时，该方军旗的位置会被亮出。
"""

#: Heading of the enumerated legal-move block. Named because both the stub
#: completer and the probe leak test locate the block by it.
LEGAL_BLOCK_MARKER = "你的全部合法走法"
#: Headings of the two derived blocks, for the same reason.
FLAG_BLOCK_MARKER = "对方军旗只可能在"
BELIEF_BLOCK_MARKER = "从战斗结果可以推出的对方军衔范围"

#: Empty-camp glyph. Deliberately *not* the CLI's ``营``: that is also the rank
#: symbol for a major, so the CLI renders both an empty camp and a major as the
#: same character and relies on the brackets to disambiguate. A human reading a
#: coloured terminal copes; there is no reason to hand a model that collision.
CAMP_GLYPH = "※"

_LEGEND = (
    "图例：<>=你的棋子  [?]=敌方棋子（军衔未知）  [旗]=敌方军旗（已亮出）  "
    f"{CAMP_GLYPH}=空行营（安全格）  ◎=空铁路  ·=空地"
)


def _grid(obs: Observation) -> str:
    own = {piece.square: piece.kind for piece in obs.own}
    enemy = {square.square: square.revealed for square in obs.enemy}
    lines = ["      A    B    C    D    E"]
    for row in range(12):
        cells: list[str] = []
        for column in range(5):
            square = (row, column)
            if square in own:
                cells.append(f"<{SYMBOLS[own[square]]}>")
            elif square in enemy:
                cells.append("[旗]" if enemy[square] else "[?]")
            elif square in CAMPS:
                cells.append(f" {CAMP_GLYPH} ")
            elif square in RAILWAYS:
                cells.append(" ◎ ")
            else:
                cells.append(" · ")
        suffix = (
            "  大本营"
            if any((row, column) in HEADQUARTERS for column in range(5))
            else ""
        )
        lines.append(f"{row + 1:>2}  " + "  ".join(cells) + suffix)
    return "\n".join(lines)


def _battle_log(obs: Observation) -> str:
    if not obs.log:
        return "（还没有走过棋）"
    lines: list[str] = []
    for item in obs.log:
        who = "你" if item.by_me else "对方"
        move = f"{format_position(item.src)}-{format_position(item.dst)}"
        if not item.had_battle:
            lines.append(f"  第{item.ply + 1}手 {who}走 {move}")
            continue
        if item.outcome is None:
            result = "结果未知"
        elif item.outcome > 0:
            result = "攻方获胜"
        elif item.outcome < 0:
            result = "守方获胜"
        else:
            result = "同归于尽"
        mine = (
            f"，你参战的是{SYMBOLS[item.own_kind]}"
            if item.own_kind is not None
            else ""
        )
        lines.append(f"  第{item.ply + 1}手 {who}攻 {move}：{result}{mine}")
    return "\n".join(lines)


MOVE_INSTRUCTION = (
    "请选择一步棋。可以先简短分析，但**最后一行只输出走法本身**，"
    "格式为 起点-终点，例如：A10-A9"
)


def render(
    obs: Observation, scaffold: Scaffold, instruction: str = MOVE_INSTRUCTION
) -> str:
    """Turn an observation into the prompt body.

    Reads nothing but ``obs``, which is what makes the leak guarantee hold: if
    a fact is not on the observation it cannot reach the model. ``instruction``
    is the trailing ask -- the diagnostic battery swaps in its own question and
    otherwise reuses this renderer unchanged.
    """
    blocks: list[str] = []
    if scaffold.rules:
        blocks.append(RULES)

    mine = SIDE_NAMES[obs.me]
    theirs = SIDE_NAMES[obs.me.other]
    blocks.append(
        f"你执{mine}，你的棋子起始于第 {SIDE_ROWS[obs.me]} 行；"
        f"对手执{theirs}，起始于第 {SIDE_ROWS[obs.me.other]} 行。\n"
        f"当前是第 {obs.ply + 1} 手，轮到你走。"
    )
    blocks.append(_LEGEND + "\n" + _grid(obs))

    if scaffold.piece_list:
        listed = "  ".join(
            f"{format_position(piece.square)}={SYMBOLS[piece.kind]}"
            for piece in obs.own
        )
        blocks.append(f"你的棋子（{len(obs.own)} 枚）：\n  {listed}")

    blocks.append("行棋记录：\n" + _battle_log(obs))

    if scaffold.flag_candidates:
        candidates = "、".join(format_position(s) for s in obs.enemy_flag_squares)
        blocks.append(
            f"推断：{FLAG_BLOCK_MARKER}这些大本营格之一 —— "
            f"{candidates or '（无）'}。\n"
            "（依据：大本营的棋子永不移动，军旗一定在大本营，"
            "所以还站着己方棋子的大本营格就是全部可能位置。）"
        )

    if scaffold.belief and obs.belief:
        lines = [
            f"  {format_position(square)} 只可能是："
            + "、".join(SYMBOLS[kind] for kind in kinds)
            for square, kinds in obs.belief
        ]
        blocks.append(f"{BELIEF_BLOCK_MARKER}：\n" + "\n".join(lines))

    if scaffold.legal_moves:
        listed = "  ".join(str(move) for move in obs.legal)
        blocks.append(
            f"{LEGAL_BLOCK_MARKER}（共 {len(obs.legal)} 个）：\n  {listed}"
        )

    blocks.append(instruction)
    return "\n\n".join(blocks)
