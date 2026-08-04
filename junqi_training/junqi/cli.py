from __future__ import annotations

import argparse
import random
from pathlib import Path

from .board import CAMPS, HEADQUARTERS, RAILWAYS
from .bot import BotWeights
from .deployment import (
    load_deployment,
    random_deployment,
    save_deployment,
    strategic_deployment,
    swap_pieces,
    validate_deployment,
)
from .game import Game
from .search_bot import SearchBot
from .opponents import OracleSearchBot
from .web_export import fingerprint
from .training import DEFAULT_MODEL
from .types import (
    Owner,
    PieceKind,
    SYMBOLS,
    format_position,
    parse_move,
    parse_position,
)


def render(game: Game) -> str:
    lines = ["      A    B    C    D    E"]
    for row in range(12):
        cells: list[str] = []
        for column in range(5):
            position = (row, column)
            piece = game.board.get(position)
            if piece is None:
                marker = "营" if position in CAMPS else "◎" if position in RAILWAYS else "·"
                cells.append(f" {marker} ")
            elif piece.owner == Owner.BOT:
                cells.append(
                    "[旗]"
                    if piece.kind == PieceKind.FLAG and piece.revealed
                    else "[?]"
                )
            else:
                symbol = SYMBOLS[piece.kind]
                brackets = ("<", ">") if piece.owner == Owner.HUMAN else ("[", "]")
                cells.append(f"{brackets[0]}{symbol}{brackets[1]}")
        hq = (
            "  大本营"
            if any((row, column) in HEADQUARTERS for column in range(5))
            else ""
        )
        lines.append(f"{row + 1:>2}  " + "  ".join(cells) + hq)
    return "\n".join(lines)


def arrange_player(game: Game, rng: random.Random) -> bool:
    print("\n=== 排兵布阵 ===")
    print("输入 swap A12 B12 交换两枚棋子；random 重新随机；done 确认；quit 退出。")
    print("save my-opening.json 保存当前阵型，下次用 --my-deployment 直接载入。")
    print("限制：军旗在大本营，地雷在最后两排，炸弹不能在最前排。\n")
    while True:
        print(render(game))
        raw = input("\n布阵 > ").strip()
        command = raw.lower()
        if command in {"done", "d", "开始"}:
            errors = validate_deployment(game.board, Owner.HUMAN)
            if errors:
                print("当前阵型不合法：" + "；".join(errors))
                continue
            print("布阵完成，对局开始！\n")
            return True
        if command in {"random", "r", "随机"}:
            game.board = {
                position: piece
                for position, piece in game.board.items()
                if piece.owner != Owner.HUMAN
            }
            game.board.update(random_deployment(Owner.HUMAN, rng))
            print("已重新生成合法随机阵型。\n")
            continue
        if command in {"quit", "q", "exit"}:
            print("已退出。")
            return False
        parts = raw.replace(",", " ").split()
        if parts and parts[0].lower() in {"save", "保存"}:
            if len(parts) != 2:
                print("无法保存：格式应为：save my-opening.json\n")
                continue
            try:
                save_deployment(game.board, Path(parts[1]), Owner.HUMAN)
            except (ValueError, OSError) as error:
                print(f"无法保存：{error}\n")
                continue
            print(
                f"已保存到 {parts[1]}。下次用 "
                f"--my-deployment {parts[1]} 直接载入。\n"
            )
            continue
        try:
            if len(parts) != 3 or parts[0].lower() not in {"swap", "s", "换"}:
                raise ValueError("格式应为：swap A12 B12")
            left, right = parse_position(parts[1]), parse_position(parts[2])
            swap_pieces(game.board, Owner.HUMAN, left, right)
            print(f"已交换 {parts[1].upper()} 与 {parts[2].upper()}。\n")
        except ValueError as error:
            print(f"无法调整：{error}\n")



#: Symbols for a disclosed board, matching the browser replay so one parser
#: reads both.
def _grid(board, rows: range, owner: Owner) -> list[str]:
    from .board import CAMPS

    lines = ["     A  B  C  D  E"]
    for row in rows:
        cells = []
        for column in range(5):
            square = (row, column)
            piece = board.get(square)
            if square in CAMPS:
                cells.append(" ·")
            elif piece is None:
                cells.append("  ")
            else:
                cells.append(f" {SYMBOLS[piece.kind]}")
        lines.append(f"{row + 1:>3}{''.join(cells)}")
    return lines


def format_replay(game: Game, opening: dict, label: str, weights: BotWeights) -> str:
    """A pasteable trajectory in the browser's replay format.

    Deliberately byte-compatible with `web/lib/replay.ts` output so the same
    analysis reads either. Both openings are disclosed only here, at the end --
    the point of a replay is to be readable *after* the game.
    """
    if game.winner == Owner.HUMAN:
        result = f"human wins at ply {game.move_count}"
    elif game.winner == Owner.BOT:
        result = f"bot wins at ply {game.move_count}"
    else:
        result = f"draw at ply {game.move_count}"

    out = [
        "JQ/60 replay",
        f"result:     {result}",
        f"difficulty: {label}",
        "engine:     python cli",
        f"weights:    {fingerprint(weights)}",
        "",
        "bot opening (rows 1-6, disclosed now the game is over):",
        *_grid(opening, range(0, 6), Owner.BOT),
        "",
        "your opening (rows 7-12):",
        *_grid(opening, range(6, 12), Owner.HUMAN),
        "",
        "moves:",
    ]
    for index, record in enumerate(game.records, start=1):
        side = "Y" if record.attacker.owner == Owner.HUMAN else "B"
        notes = []
        if record.defender is not None:
            if record.defender.kind == PieceKind.FLAG:
                notes.append("军旗被夺")
            elif record.outcome > 0:
                notes.append("攻方胜")
            elif record.outcome < 0:
                notes.append("守方胜")
            else:
                notes.append("同归于尽")
        move = record.move
        text = (
            f"{index:>3} {side}  {move}  "
            f"{format_position(move.src)} → {format_position(move.dst)}"
        )
        if notes:
            text += "  · " + "  · ".join(notes)
        out.append(text)
    return "\n".join(out) + "\n"


def play(
    seed: int | None,
    model_path: Path,
    auto_deploy: bool = False,
    search_samples: int = 6,
    deployment_model: Path | None = None,
    search_replies: int = 4,
    opponent: str = "search",
    replay_path: Path | None = None,
    player_deployment: Path | None = None,
) -> None:
    weights = BotWeights.load(model_path) if model_path.exists() else BotWeights()
    if opponent == "search":
        bot = SearchBot(
            weights, seed=seed, samples=search_samples, reply_width=search_replies
        )
    else:
        # The cheating opponent. It reads every one of your hidden ranks; the
        # point is to find out how much that is actually worth against a person.
        bot = OracleSearchBot(
            weights,
            seed=seed,
            samples=search_samples,
            reply_width=search_replies,
            gamma=0.5 if opponent == "oracle-half" else 1.0,
            gamma_beam=1.0 if opponent == "oracle-perfect" else 0.0,
        )
    rng = random.Random(seed)
    bot_deployment = (
        load_deployment(deployment_model)
        if deployment_model is not None and deployment_model.exists()
        else strategic_deployment(Owner.BOT, rng)
    )
    game = Game.new(seed=seed, bot_deployment=bot_deployment)
    if player_deployment is not None:
        # Replace the random human layout with the saved one, *before*
        # arranging, so it can still be tweaked and re-saved. A bad file is
        # fatal rather than silently ignored: playing a random opening you
        # believed was yours is worse than being told to fix the path.
        loaded = load_deployment(player_deployment, Owner.HUMAN)
        game.board = {
            position: piece
            for position, piece in game.board.items()
            if piece.owner != Owner.HUMAN
        }
        game.board.update(loaded)
        print(f"已载入你的阵型：{player_deployment}")
    if opponent != "search":
        print(
            f"*** {opponent.upper()}: this opponent SEES YOUR HIDDEN RANKS. "
            f"It is a measuring device, not the shipped bot. ***"
        )
    print("军棋单人版：你执南方（<棋>），Bot 执北方（[?]）。")
    if not auto_deploy and not arrange_player(game, rng):
        return
    # Snapshot AFTER arranging, so the replay discloses the board actually
    # played rather than the pre-swap one.
    opening = dict(game.board)
    print("输入 A10-A9 走棋；输入 moves 查看合法走法；输入 quit 退出。\n")

    while not game.over:
        print(render(game))
        if game.turn == Owner.HUMAN:
            raw = input("\n你的回合 > ").strip()
            if raw.lower() in {"quit", "q", "exit"}:
                print("已退出本局。")
                return
            if raw.lower() in {"moves", "m"}:
                print("合法走法：" + "  ".join(map(str, game.legal_moves())))
                continue
            try:
                message = game.apply(parse_move(raw))
            except ValueError as error:
                print(f"无法走棋：{error}")
                continue
            print(message)
        else:
            move = bot.choose_move(game)
            print(f"\nBot 走：{move}")
            print(game.apply(move))
        print()

    print(render(game))
    if game.winner == Owner.HUMAN:
        print("\n你赢了！")
    elif game.winner == Owner.BOT:
        print("\nBot 获胜。")
    else:
        print("\n和棋。")

    if replay_path is not None:
        replay_path.write_text(
            format_replay(game, opening, opponent, weights),
            encoding="utf-8",
            newline="\n",
        )
        print(f"\n棋谱已写入 {replay_path}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="军棋单人对局")
    parser.add_argument("--seed", type=int, default=None, help="固定棋局随机种子")
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL, help="Bot 模型路径")
    parser.add_argument(
        "--auto-deploy",
        action="store_true",
        help="跳过手工布阵，直接使用随机合法阵型",
    )
    parser.add_argument(
        "--search-samples",
        type=int,
        default=6,
        help="Bot 每个候选走法的暗子采样数；越大越强但思考更久",
    )
    parser.add_argument(
        "--search-replies",
        type=int,
        default=4,
        help="Bot 每次推演考虑的对手最强应手数；越大越谨慎但思考更久",
    )
    parser.add_argument(
        "--deployment-model",
        type=Path,
        default=None,
        help="固定的 Bot 初始布局模型路径；默认每局重新生成战术布阵",
    )
    parser.add_argument(
        "--opponent",
        choices=("search", "oracle", "oracle-half", "oracle-perfect"),
        default="search",
        help=(
            "对手类型。search 是正常出厂 Bot；"
            "oracle **会看穿你的全部暗子**（同样的权重，只是不再采样猜测），"
            "但候选走法仍由不作弊的启发式挑选；"
            "oracle-perfect 连挑candidate也用真实军衔，是真正的上界；"
            "oracle-half 只看穿一半。这三个是测量工具，不是会发布的对手"
        ),
    )
    parser.add_argument(
        "--replay",
        type=Path,
        default=None,
        help="终局后把棋谱写到这个文件，格式与网页版一致，可直接贴出来分析",
    )
    parser.add_argument(
        "--my-deployment",
        type=Path,
        default=None,
        help=(
            "载入你自己的初始阵型（布阵阶段用 save my-opening.json 保存的文件）。"
            "载入后仍可继续 swap 调整再 done；文件不合法会直接报错退出，"
            "不会悄悄换成随机阵型"
        ),
    )
    return parser


def main() -> None:
    arguments = build_parser().parse_args()
    if arguments.search_samples < 1:
        raise SystemExit("search-samples 至少为 1")
    if arguments.search_replies < 1:
        raise SystemExit("search-replies 至少为 1")
    if arguments.my_deployment is not None and not arguments.my_deployment.exists():
        raise SystemExit(f"--my-deployment {arguments.my_deployment} 不存在")
    if arguments.my_deployment is not None and arguments.auto_deploy:
        raise SystemExit("--my-deployment 与 --auto-deploy 冲突：前者就是你的阵型")
    play(
        arguments.seed,
        arguments.model,
        arguments.auto_deploy,
        arguments.search_samples,
        arguments.deployment_model,
        arguments.search_replies,
        arguments.opponent,
        arguments.replay,
        arguments.my_deployment,
    )


if __name__ == "__main__":
    main()
