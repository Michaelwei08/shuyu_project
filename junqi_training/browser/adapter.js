/**
 * Point our engine at junqi.app's board.
 *
 * Loaded into the page alongside dist-engine/junqi-engine.js. Reads their DOM,
 * hands the visible state to our information-set search, and clicks the move.
 *
 * Board contract observed on their page:
 *   60 elements matching `.Square`, on a 12x5 grid of centre points.
 *   class contains "Blue"        -> our piece, rank is the element text
 *   class contains "Red"         -> their piece, rank not shown (concealed)
 *   class contains "Transparent" -> empty square
 *
 * We sit at the bottom of their board and our engine plays the top ("bot"),
 * so every square is rotated 180 degrees. The board is symmetric under that
 * rotation -- camps, headquarters and railways all map onto themselves -- so
 * the transform is exact rather than an approximation.
 *
 * Known handicap: their DOM does not expose structured battle results, so the
 * engine starts each turn with an empty belief instead of accumulating the
 * rank deductions it normally makes from past fights.
 */
(function () {
  const LABEL_TO_KIND = {
    "司令": "COMMANDER", "军长": "GENERAL", "师长": "MAJOR_GENERAL",
    "旅长": "BRIGADIER", "团长": "COLONEL", "营长": "MAJOR",
    "连长": "CAPTAIN", "排长": "LIEUTENANT", "工兵": "ENGINEER",
    "地雷": "MINE", "炸弹": "BOMB", "军旗": "FLAG",
  };
  const ROWS = 12;
  const COLS = 5;

  function axis(values, count) {
    const sorted = [...values].sort((a, b) => a - b);
    const groups = [];
    for (const value of sorted) {
      const last = groups[groups.length - 1];
      if (last && Math.abs(last[0] - value) < 25) last.push(value);
      else groups.push([value]);
    }
    if (groups.length !== count) {
      throw new Error(`expected ${count} lines, found ${groups.length}`);
    }
    return groups.map((g) => g.reduce((a, b) => a + b, 0) / g.length);
  }

  /** Every square, indexed in OUR coordinate frame (0 = A1, north-west). */
  function scan() {
    const nodes = [...document.querySelectorAll(".Square")];
    if (nodes.length !== ROWS * COLS) {
      throw new Error(`expected 60 squares, found ${nodes.length}`);
    }
    const boxes = nodes.map((el) => {
      const r = el.getBoundingClientRect();
      return { el, x: r.x + r.width / 2, y: r.y + r.height / 2 };
    });
    const xs = axis(boxes.map((b) => b.x), COLS);
    const ys = axis(boxes.map((b) => b.y), ROWS);
    const nearest = (v, lines) =>
      lines.reduce((best, l, i) =>
        Math.abs(l - v) < Math.abs(lines[best] - v) ? i : best, 0);

    const cells = new Array(ROWS * COLS).fill(null);
    for (const box of boxes) {
      const screenRow = nearest(box.y, ys);
      const screenCol = nearest(box.x, xs);
      // 180 degree rotation: we play the bottom, the engine plays the top.
      const index = (ROWS - 1 - screenRow) * COLS + (COLS - 1 - screenCol);
      const cls = (box.el.className || "").toString();
      cells[index] = {
        el: box.el,
        side: cls.includes("Blue") ? "bot" : cls.includes("Red") ? "human" : null,
        label: box.el.textContent.trim(),
      };
    }
    return cells;
  }

  const coordinate = (index) =>
    `${"ABCDE"[index % COLS]}${Math.floor(index / COLS) + 1}`;

  function readBoard() {
    const board = {};
    scan().forEach((cell, index) => {
      if (!cell || !cell.side) return;
      const entry = { owner: cell.side };
      if (cell.side === "bot") {
        const kind = LABEL_TO_KIND[cell.label];
        if (!kind) throw new Error(`unknown rank label: "${cell.label}"`);
        entry.kind = kind;
      }
      board[coordinate(index)] = entry;
    });
    return board;
  }

  const statusText = () => (document.body.innerText.match(/现在是.*?回合[！!]?/) || [])[0] || "";
  const ourTurn = () => statusText().includes("蓝方");
  const finished = () => /获胜|胜利|结束|失败|赢/.test(document.body.innerText.slice(0, 400));

  function squareAt(coord) {
    const col = "ABCDE".indexOf(coord[0]);
    const row = Number(coord.slice(1)) - 1;
    return scan()[row * COLS + col].el;
  }

  const trajectory = [];

  async function step(difficulty) {
    if (!ourTurn()) return { acted: false, reason: "not our turn" };
    const board = readBoard();
    const move = window.JunqiEngine.chooseMove(board, difficulty);
    squareAt(move.from).click();
    await new Promise((r) => setTimeout(r, 220));
    squareAt(move.to).click();
    trajectory.push(`${trajectory.length + 1} ${move.from}-${move.to}`);
    return { acted: true, move };
  }

  window.__junqiAdapter = {
    readBoard,
    ourTurn,
    finished,
    statusText,
    trajectory,
    step,
    /** Play until it is no longer our move to make. */
    async run(maxPlies = 200, difficulty = "focused") {
      const log = [];
      for (let i = 0; i < maxPlies; i += 1) {
        if (finished()) { log.push("game over: " + statusText()); break; }
        if (!ourTurn()) { await new Promise((r) => setTimeout(r, 700)); continue; }
        try {
          const result = await step(difficulty);
          if (result.acted) log.push(`${result.move.from}-${result.move.to}`);
        } catch (error) {
          log.push("ERROR " + error.message);
          break;
        }
        await new Promise((r) => setTimeout(r, 900));
      }
      return { plies: trajectory.length, tail: log.slice(-8), status: statusText() };
    },
  };
  return "adapter ready";
})();
