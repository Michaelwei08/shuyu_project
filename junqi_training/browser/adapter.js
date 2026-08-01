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

  /**
   * Every square, indexed in OUR coordinate frame (0 = A1, north-west).
   *
   * Rank into fixed 12x5 chunks rather than clustering coordinates by
   * tolerance: the page relayouts mid-game (column centres shifted from
   * 591..1117 to 608..1099 during the first run) and a transient frame briefly
   * produced a sixth column, which killed the loop.
   */
  function scan() {
    const nodes = [...document.querySelectorAll(".Square")];
    if (nodes.length !== ROWS * COLS) {
      throw new Error(`expected ${ROWS * COLS} squares, found ${nodes.length}`);
    }
    const boxes = nodes
      .map((el) => {
        const r = el.getBoundingClientRect();
        return { el, x: r.x + r.width / 2, y: r.y + r.height / 2 };
      })
      .sort((a, b) => a.y - b.y || a.x - b.x);

    const cells = new Array(ROWS * COLS).fill(null);
    for (let screenRow = 0; screenRow < ROWS; screenRow += 1) {
      const row = boxes
        .slice(screenRow * COLS, screenRow * COLS + COLS)
        .sort((a, b) => a.x - b.x);
      row.forEach((box, screenCol) => {
        // 180 degree rotation: we play the bottom, the engine plays the top.
        const index = (ROWS - 1 - screenRow) * COLS + (COLS - 1 - screenCol);
        const cls = (box.el.className || "").toString();
        cells[index] = {
          el: box.el,
          side: cls.includes("Blue") ? "bot" : cls.includes("Red") ? "human" : null,
          label: box.el.textContent.trim(),
        };
      });
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

  const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

  /** Our own ranks, keyed by square, in the engine's frame. */
  function ourLayout() {
    const map = {};
    scan().forEach((cell, index) => {
      if (cell && cell.side === "bot") map[coordinate(index)] = LABEL_TO_KIND[cell.label];
    });
    return map;
  }

  /**
   * Replace their default opening with one from our generator.
   *
   * Their deployment UI swaps two of your own pieces per click-pair, same as
   * ours, so any legal target layout is reachable by a permutation sort. This
   * matters: the generator screens the flag headquarters with all three mines,
   * which is the defence that stretched a real game from ply 37 to 57.
   */
  async function installDeployment() {
    const target = window.JunqiEngine.suggestDeployment();
    // Their UI validates every swap, so an intermediate that puts a mine off
    // the rear rows, a bomb on the front row or the flag outside a
    // headquarters is simply refused -- 15 of 23 were, in the first attempt.
    // Settle the constrained ranks first: once all three mines sit on their
    // targets no later swap can displace one, and the same for flag and bombs.
    // Headquarters squares go first, both of them. A piece parked in a
    // headquarters can never move again, so a half-applied layout that strands
    // something valuable there is worse than either the source or the target:
    // a real game ended up with a 军长 (general) frozen in one.
    const priority = { FLAG: 0, MINE: 1, BOMB: 2 };
    const HEADQUARTERS_SQUARES = ["B1", "D1"];
    const rank = (square) =>
      HEADQUARTERS_SQUARES.includes(square) ? -1 : (priority[target[square]] ?? 3);
    const squares = Object.keys(target).sort((a, b) => rank(a) - rank(b));
    let swaps = 0;
    let refused = 0;
    for (const square of squares) {
      const now = ourLayout();
      if (now[square] === target[square]) continue;
      const donor = squares.find(
        (s) => s !== square && now[s] === target[square] && now[s] !== target[s],
      );
      if (!donor) { refused += 1; continue; }
      squareAt(square).click();
      await sleep(140);
      squareAt(donor).click();
      await sleep(200);
      if (ourLayout()[square] === target[square]) swaps += 1;
      else refused += 1;
    }
    // Last resort: whatever ended up frozen in the spare headquarters must be
    // cheap, even if the rest of the layout could not be reached.
    const EXPENDABLE = ["LIEUTENANT", "CAPTAIN", "ENGINEER", "MAJOR"];
    for (const hq of HEADQUARTERS_SQUARES) {
      const now = ourLayout();
      if (now[hq] === "FLAG" || EXPENDABLE.includes(now[hq])) continue;
      const donor = Object.keys(now).find(
        (s) => !HEADQUARTERS_SQUARES.includes(s) && EXPENDABLE.includes(now[s]),
      );
      if (!donor) continue;
      squareAt(hq).click();
      await sleep(140);
      squareAt(donor).click();
      await sleep(200);
      if (EXPENDABLE.includes(ourLayout()[hq])) swaps += 1;
    }

    const final = ourLayout();
    const wrong = squares.filter((s) => final[s] !== target[s]);
    return {
      swaps,
      refused,
      mismatched: wrong.length,
      headquarters: HEADQUARTERS_SQUARES.map((s) => `${s}=${final[s]}`),
    };
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

  const recent = [];

  /**
   * What past battles proved about each enemy square.
   *
   * Their DOM has no structured battle log, so this is rebuilt by diffing the
   * board around each of our own attacks -- we know our attacker's rank, and
   * the diff tells us whether we won, traded, or died. Without it the engine
   * re-attacks a square that just killed a stronger piece, which is exactly
   * the "connect died, now try the major" blunder: any rank that beat our
   * captain either beats or ties the major too.
   */
  const beliefs = {};
  let lastSeen = null;

  const redSquares = (board) =>
    Object.keys(board).filter((sq) => board[sq].owner === "human").sort();

  /** Follow a belief when its piece walks, drop it when the square resolves. */
  function reconcileOpponentMoves(board) {
    if (!lastSeen) { lastSeen = redSquares(board); return; }
    const now = redSquares(board);
    const vanished = lastSeen.filter((s) => !now.includes(s));
    const appeared = now.filter((s) => !lastSeen.includes(s));
    if (vanished.length === 1 && appeared.length === 1) {
      if (beliefs[vanished[0]]) beliefs[appeared[0]] = beliefs[vanished[0]];
    }
    for (const square of vanished) delete beliefs[square];
    lastSeen = now;
  }

  /**
   * Wait for their board to actually finish moving.
   *
   * Reading a fixed 320ms after the click was not enough: the source square
   * had cleared but the destination still showed the old occupant, so every
   * single attack was classified "lost". Over one game that produced 21
   * battles with zero wins while material was in fact an even exchange, and
   * the false constraints it wrote into `beliefs` made the bot re-attack
   * squares it had already taken. Require two consecutive identical reads that
   * differ from the starting position.
   */
  async function settle(snapshot, timeoutMs = 3000) {
    const deadline = Date.now() + timeoutMs;
    let previous = null;
    while (Date.now() < deadline) {
      await sleep(120);
      const now = readBoard();
      const serialised = JSON.stringify(now);
      if (serialised !== snapshot && serialised === previous) return now;
      previous = serialised;
    }
    return readBoard();
  }

  function learnFromAttack(attackerKind, target, before, after) {
    const wasEnemy = before[target] && before[target].owner === "human";
    if (!wasEnemy) return null;
    const nowOurs = after[target] && after[target].owner === "bot";
    if (nowOurs) { delete beliefs[target]; return "won"; }
    if (!after[target]) { delete beliefs[target]; return "traded"; }
    // Their piece held the square, so it beats our attacker. Intersect with
    // anything we already knew about it.
    const survivors = window.JunqiEngine.survivorsOf(attackerKind);
    beliefs[target] = beliefs[target]
      ? beliefs[target].filter((k) => survivors.includes(k))
      : survivors;
    return "lost";
  }

  /**
   * The engine is handed a fresh state each turn, because their DOM exposes no
   * structured battle log to rebuild belief from. With no history it cannot
   * see that it is repeating itself, and it oscillated D8-D9-D8-D9 for six
   * plies in a real game. Reject a move that simply undoes the last one, and
   * let the search's noise produce an alternative.
   */
  function undoesLastMove(move) {
    const last = recent[recent.length - 1];
    return Boolean(last && last.from === move.to && last.to === move.from);
  }

  async function step(difficulty) {
    if (!ourTurn()) return { acted: false, reason: "not our turn" };
    const before = readBoard();
    reconcileOpponentMoves(before);
    let move = window.JunqiEngine.chooseMove(before, difficulty, beliefs);
    for (let retry = 0; retry < 4 && undoesLastMove(move); retry += 1) {
      move = window.JunqiEngine.chooseMove(before, difficulty, beliefs);
    }
    const attackerKind = before[move.from] && before[move.from].kind;
    const snapshot = JSON.stringify(before);
    squareAt(move.from).click();
    await sleep(200);
    squareAt(move.to).click();
    const after = await settle(snapshot);
    // A click the page ignores would otherwise be re-issued forever.
    const landed = JSON.stringify(after) !== snapshot;
    let outcome = null;
    if (landed) {
      outcome = learnFromAttack(attackerKind, move.to, before, after);
      lastSeen = redSquares(after);
      recent.push(move);
      if (recent.length > 8) recent.shift();
      trajectory.push(
        `${trajectory.length + 1} ${move.from}-${move.to}` +
        (outcome ? ` (${outcome})` : ""),
      );
    }
    return { acted: landed, move, landed, outcome, beliefs: Object.keys(beliefs).length };
  }

  window.__junqiAdapter = {
    readBoard,
    ourLayout,
    installDeployment,
    ourTurn,
    finished,
    statusText,
    trajectory,
    beliefs,
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
