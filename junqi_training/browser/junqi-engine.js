var JunqiEngine = (function(exports) {
	Object.defineProperty(exports, Symbol.toStringTag, { value: "Module" });
	//#region lib/board.ts
	var ROWS = 12;
	var COLS = 5;
	var indexOf = (row, col) => row * 5 + col;
	var rowOf = (index) => Math.floor(index / 5);
	var colOf = (index) => index % 5;
	var CAMPS = new Set([
		indexOf(2, 1),
		indexOf(2, 3),
		indexOf(3, 2),
		indexOf(4, 1),
		indexOf(4, 3),
		indexOf(7, 1),
		indexOf(7, 3),
		indexOf(8, 2),
		indexOf(9, 1),
		indexOf(9, 3)
	]);
	var HEADQUARTERS = new Set([
		indexOf(0, 1),
		indexOf(0, 3),
		indexOf(11, 1),
		indexOf(11, 3)
	]);
	var RAILS = new Set([
		...[
			1,
			5,
			6,
			10
		].flatMap((row) => Array.from({ length: 5 }, (_, col) => indexOf(row, col))),
		...Array.from({ length: 10 }, (_, offset) => offset + 1).flatMap((row) => [indexOf(row, 0), indexOf(row, 4)]),
		indexOf(5, 2),
		indexOf(6, 2)
	]);
	var RIVER_FILES = new Set([
		0,
		2,
		4
	]);
	var ORTHOGONAL_STEPS = [
		[1, 0],
		[-1, 0],
		[0, 1],
		[0, -1]
	];
	var UNREACHABLE = 60;
	function inBounds(row, col) {
		return row >= 0 && row < 12 && col >= 0 && col < 5;
	}
	function crossesBlockedRiver(fromRow, toRow, col) {
		return Math.min(fromRow, toRow) === 5 && Math.max(fromRow, toRow) === 6 && !RIVER_FILES.has(col);
	}
	/** Static board topology. Occupancy is applied later in `legalMoves`. */
	var ROAD_NEIGHBORS = Array.from({ length: 60 }, (_, index) => {
		const row = rowOf(index);
		const col = colOf(index);
		const result = [];
		for (const [dr, dc] of ORTHOGONAL_STEPS) {
			const nextRow = row + dr;
			const nextCol = col + dc;
			if (inBounds(nextRow, nextCol) && !crossesBlockedRiver(row, nextRow, col)) result.push(indexOf(nextRow, nextCol));
		}
		if (CAMPS.has(index)) {
			for (const [dr, dc] of [
				[1, 1],
				[1, -1],
				[-1, 1],
				[-1, -1]
			]) if (inBounds(row + dr, col + dc)) result.push(indexOf(row + dr, col + dc));
		}
		for (const camp of CAMPS) if (Math.abs(row - rowOf(camp)) === 1 && Math.abs(col - colOf(camp)) === 1) result.push(camp);
		return [...new Set(result)];
	});
	var RAIL_NEIGHBORS = Array.from({ length: 60 }, (_, index) => {
		if (!RAILS.has(index)) return [];
		const row = rowOf(index);
		const col = colOf(index);
		return ORTHOGONAL_STEPS.flatMap(([dr, dc]) => {
			const nextRow = row + dr;
			const nextCol = col + dc;
			if (!inBounds(nextRow, nextCol) || crossesBlockedRiver(row, nextRow, col)) return [];
			const next = indexOf(nextRow, nextCol);
			return RAILS.has(next) ? [next] : [];
		});
	});
	var RAIL_RAYS = Array.from({ length: 60 }, (_, index) => {
		if (!RAILS.has(index)) return [];
		const row = rowOf(index);
		const col = colOf(index);
		return ORTHOGONAL_STEPS.map(([dr, dc]) => {
			const ray = [];
			let nextRow = row + dr;
			let nextCol = col + dc;
			while (inBounds(nextRow, nextCol)) {
				const next = indexOf(nextRow, nextCol);
				if (!RAILS.has(next) || crossesBlockedRiver(nextRow - dr, nextRow, nextCol)) break;
				ray.push(next);
				nextRow += dr;
				nextCol += dc;
			}
			return ray;
		}).filter((ray) => ray.length);
	});
	/**
	* Moves -- not squares -- between every pair, on an empty board.
	*
	* Manhattan distance is badly wrong here: a piece on a railway crosses the
	* whole rank in one move, so E2 is two moves from B1 while being Manhattan 4
	* away. Blockers only slow a piece down, so the empty board is a lower bound --
	* the conservative direction for judging danger. Nothing is more than 5 moves
	* from anything else.
	*/
	var MOVE_DISTANCE = (() => {
		const empty = Array(60).fill(null);
		return Array.from({ length: 60 }, (_, source) => {
			const distances = Array(60).fill(UNREACHABLE);
			distances[source] = 0;
			const queue = [source];
			for (let head = 0; head < queue.length; head += 1) {
				const current = queue[head];
				const step = distances[current] + 1;
				const reachable = new Set([...ROAD_NEIGHBORS[current], ...railDestinations(current, "CAPTAIN", empty)]);
				for (const next of reachable) if (distances[next] === UNREACHABLE) {
					distances[next] = step;
					queue.push(next);
				}
			}
			return distances;
		});
	})();
	function moveDistance(source, target) {
		return MOVE_DISTANCE[source][target];
	}
	/**
	* True when only an engineer could make this move.
	*
	* Engineers alone turn corners on the railway, so a move reachable by the
	* engineer's rail BFS but not by a step or a straight slide announces the
	* piece's rank to anyone watching. Spending that disguise for nothing is how
	* an engineer gets picked off.
	*/
	function revealsEngineer(board, from, to) {
		if (ROAD_NEIGHBORS[from].includes(to)) return false;
		if (railDestinations(from, "CAPTAIN", board).includes(to)) return false;
		return railDestinations(from, "ENGINEER", board).includes(to);
	}
	function coordinate(index) {
		return `${String.fromCharCode(65 + colOf(index))}${rowOf(index) + 1}`;
	}
	function deploymentSquares(owner) {
		return (owner === "bot" ? Array.from({ length: 6 }, (_, i) => i) : Array.from({ length: 6 }, (_, i) => i + 6)).flatMap((row) => Array.from({ length: 5 }, (_, col) => indexOf(row, col))).filter((index) => !CAMPS.has(index));
	}
	function ownHeadquarters(owner) {
		return owner === "bot" ? [indexOf(0, 1), indexOf(0, 3)] : [indexOf(11, 1), indexOf(11, 3)];
	}
	function rearRows(owner) {
		return owner === "bot" ? BOT_REAR_ROWS : HUMAN_REAR_ROWS;
	}
	function frontRow(owner) {
		return owner === "bot" ? 5 : 6;
	}
	var BOT_REAR_ROWS = new Set([0, 1]);
	var HUMAN_REAR_ROWS = new Set([10, 11]);
	function railDestinations(from, kind, board) {
		if (!RAILS.has(from)) return [];
		if (kind === "ENGINEER") {
			const visited = new Set([from]);
			const queue = [from];
			const result = [];
			while (queue.length) {
				const current = queue.shift();
				for (const next of RAIL_NEIGHBORS[current]) {
					if (visited.has(next)) continue;
					visited.add(next);
					result.push(next);
					if (!board[next]) queue.push(next);
				}
			}
			return result;
		}
		const result = [];
		for (const ray of RAIL_RAYS[from]) for (const next of ray) {
			result.push(next);
			if (board[next]) break;
		}
		return result;
	}
	/**
	* Headquarters squares that may still hide `owner`'s flag.
	*
	* Deployment fills every non-camp home square, headquarters pieces can never
	* move, and a flag never leaves its headquarters -- so the flag is under one of
	* the own-side pieces still standing on a headquarters square. A revealed flag
	* collapses this to a single square. Reads occupancy only, never a rank, so
	* either side may use it against the other.
	*/
	function liveFlagSquares(board, owner, knowledge) {
		const held = ownHeadquarters(owner).filter((index) => board[index]?.owner === owner);
		const revealed = held.filter((index) => board[index]?.revealed);
		if (revealed.length) return revealed;
		const possible = knowledge ? held.filter((index) => !knowledge[index] || knowledge[index].includes("FLAG")) : held;
		return possible.length ? possible : held;
	}
	function legalMoves(board, owner) {
		const moves = [];
		board.forEach((piece, from) => {
			if (!piece || piece.owner !== owner || piece.kind === "FLAG" || piece.kind === "MINE") return;
			if (HEADQUARTERS.has(from)) return;
			const candidates = new Set([...ROAD_NEIGHBORS[from], ...railDestinations(from, piece.kind, board)]);
			for (const to of candidates) {
				const target = board[to];
				if (!target || target.owner !== owner && !CAMPS.has(to)) moves.push({
					from,
					to
				});
			}
		});
		return moves;
	}
	//#endregion
	//#region lib/types.ts
	var KINDS = [
		"FLAG",
		"COMMANDER",
		"GENERAL",
		"MAJOR_GENERAL",
		"BRIGADIER",
		"COLONEL",
		"MAJOR",
		"CAPTAIN",
		"LIEUTENANT",
		"ENGINEER",
		"MINE",
		"BOMB"
	];
	var COUNTS = {
		FLAG: 1,
		COMMANDER: 1,
		GENERAL: 1,
		MAJOR_GENERAL: 2,
		BRIGADIER: 2,
		COLONEL: 2,
		MAJOR: 2,
		CAPTAIN: 3,
		LIEUTENANT: 3,
		ENGINEER: 3,
		MINE: 3,
		BOMB: 2
	};
	var RANK = {
		FLAG: 0,
		COMMANDER: 1,
		GENERAL: 2,
		MAJOR_GENERAL: 3,
		BRIGADIER: 4,
		COLONEL: 5,
		MAJOR: 6,
		CAPTAIN: 7,
		LIEUTENANT: 8,
		ENGINEER: 9,
		MINE: 10,
		BOMB: 11
	};
	function other(owner) {
		return owner === "human" ? "bot" : "human";
	}
	//#endregion
	//#region lib/belief.ts
	function sampleHiddenState(state) {
		const sampled = {
			...state,
			board: state.board.map((piece) => piece ? { ...piece } : null),
			history: [...state.history],
			knowledge: { ...state.knowledge }
		};
		const positions = sampled.board.flatMap((piece, index) => piece?.owner === "human" && !piece.revealed ? [index] : []);
		positions.forEach((index) => {
			sampled.board[index] = null;
		});
		assignHiddenKinds(positions, state).forEach((kind, index) => {
			sampled.board[index] = {
				owner: "human",
				kind,
				revealed: false
			};
		});
		return sampled;
	}
	function assignHiddenKinds(positions, state) {
		const hiddenHeadquarters = positions.filter((position) => ownHeadquarters("human").includes(position) && (!state.knowledge[position] || state.knowledge[position].includes("FLAG")));
		const revealedFlag = state.board.findIndex((piece) => piece?.owner === "human" && piece.kind === "FLAG" && piece.revealed);
		const flagPosition = revealedFlag >= 0 ? null : hiddenHeadquarters[Math.floor(Math.random() * hiddenHeadquarters.length)];
		if (revealedFlag < 0 && flagPosition === void 0) throw new Error("Cannot construct a fair hidden-state sample: no legal flag square.");
		const allowed = (position, kind) => {
			if (kind === "FLAG") return position === flagPosition;
			if (kind === "MINE" && !rearRows("human").has(rowOf(position))) return false;
			return !state.knowledge[position] || state.knowledge[position].includes(kind);
		};
		for (let attempt = 0; attempt < 80; attempt += 1) {
			const result = /* @__PURE__ */ new Map();
			const remaining = { ...COUNTS };
			const aliveCount = positions.length + (revealedFlag >= 0 ? 1 : 0);
			const casualties = Object.values(COUNTS).reduce((sum, count) => sum + count, 0) - aliveCount;
			for (let removed = 0; removed < casualties; removed += 1) {
				const candidates = KINDS.flatMap((kind) => kind !== "FLAG" ? Array.from({ length: remaining[kind] }, () => kind) : []);
				const kind = candidates[Math.floor(Math.random() * candidates.length)];
				remaining[kind] -= 1;
			}
			if (revealedFlag >= 0) remaining.FLAG = 0;
			if (flagPosition !== null && flagPosition !== void 0) {
				result.set(flagPosition, "FLAG");
				remaining.FLAG = 0;
			}
			const available = positions.filter((position) => position !== flagPosition);
			const mineOptions = available.filter((position) => allowed(position, "MINE")).sort(() => Math.random() - .5);
			if (mineOptions.length < remaining.MINE) continue;
			for (let count = 0; count < remaining.MINE; count += 1) result.set(mineOptions[count], "MINE");
			const minePositions = new Set(mineOptions.slice(0, remaining.MINE));
			remaining.MINE = 0;
			const ordered = available.filter((position) => !minePositions.has(position)).sort((left, right) => {
				return KINDS.filter((kind) => allowed(left, kind)).length - KINDS.filter((kind) => allowed(right, kind)).length || Math.random() - .5;
			});
			let failed = false;
			for (const position of ordered) {
				const candidates = KINDS.flatMap((kind) => allowed(position, kind) ? Array.from({ length: remaining[kind] }, () => kind) : []);
				if (!candidates.length) {
					failed = true;
					break;
				}
				const kind = candidates[Math.floor(Math.random() * candidates.length)];
				result.set(position, kind);
				remaining[kind] -= 1;
			}
			if (!failed) return result;
		}
		throw new Error("Cannot construct a fair hidden-state sample from current knowledge.");
	}
	function sampleWorlds(state, count) {
		const worlds = [];
		for (let attempt = 0; attempt < count * 2 && worlds.length < count; attempt += 1) try {
			worlds.push(sampleHiddenState(state));
		} catch {}
		return worlds;
	}
	//#endregion
	//#region lib/deployment.ts
	function shuffle(items) {
		const copy = [...items];
		for (let index = copy.length - 1; index > 0; index -= 1) {
			const target = Math.floor(Math.random() * (index + 1));
			[copy[index], copy[target]] = [copy[target], copy[index]];
		}
		return copy;
	}
	function randomDeployment(owner) {
		const board = Array(60).fill(null);
		const available = deploymentSquares(owner);
		const flag = ownHeadquarters(owner)[Math.floor(Math.random() * 2)];
		board[flag] = {
			owner,
			kind: "FLAG",
			revealed: false
		};
		available.splice(available.indexOf(flag), 1);
		const place = (kind, count, allowed) => {
			shuffle(available.filter(allowed)).slice(0, count).forEach((index) => {
				board[index] = {
					owner,
					kind,
					revealed: false
				};
				available.splice(available.indexOf(index), 1);
			});
		};
		place("MINE", COUNTS.MINE, (index) => rearRows(owner).has(rowOf(index)));
		place("BOMB", COUNTS.BOMB, (index) => rowOf(index) !== frontRow(owner));
		shuffle(KINDS.flatMap((kind) => [
			"FLAG",
			"MINE",
			"BOMB"
		].includes(kind) ? [] : Array(COUNTS[kind]).fill(kind))).forEach((kind, index) => {
			board[available[index]] = {
				owner,
				kind,
				revealed: false
			};
		});
		return board;
	}
	function pick(items) {
		return items[Math.floor(Math.random() * items.length)];
	}
	/**
	* A fresh, legal, non-random-looking opening.
	*
	* A fixed opening is worth nothing once the opponent has seen it twice, and a
	* uniformly random one wastes material. This keeps the shape sensible -- mines
	* screening the flag, a cheap decoy in the unused headquarters, leaders off the
	* back rows -- while varying every game.
	*/
	function strategicDeployment(owner) {
		for (let attempt = 0; attempt < 60; attempt += 1) {
			const board = buildStrategicLayout(owner);
			if (board && validateDeployment(board, owner)) return board;
		}
		return randomDeployment(owner);
	}
	function buildStrategicLayout(owner) {
		const board = Array(60).fill(null);
		const free = new Set(deploymentSquares(owner));
		const rear = rearRows(owner);
		const front = frontRow(owner);
		const put = (index, kind) => {
			board[index] = {
				owner,
				kind,
				revealed: false
			};
			free.delete(index);
		};
		const [flagHQ, decoyHQ] = shuffle(ownHeadquarters(owner));
		put(flagHQ, "FLAG");
		put(decoyHQ, pick(["LIEUTENANT", "CAPTAIN"]));
		const rearFree = [...free].filter((index) => rear.has(rowOf(index)));
		const guards = shuffle(rearFree.filter((index) => Math.abs(rowOf(index) - rowOf(flagHQ)) + Math.abs(colOf(index) - colOf(flagHQ)) === 1)).slice(0, Math.random() < .75 ? 3 : 2);
		const mineSlots = [...guards, ...shuffle(rearFree.filter((index) => !guards.includes(index)))].slice(0, COUNTS.MINE);
		if (mineSlots.length < COUNTS.MINE) return null;
		mineSlots.forEach((index) => put(index, "MINE"));
		const midfield = shuffle([...free].filter((index) => rowOf(index) !== front && !rear.has(rowOf(index))));
		if (midfield.length < COUNTS.BOMB + 2) return null;
		midfield.slice(0, COUNTS.BOMB).forEach((index) => put(index, "BOMB"));
		const leaders = ["COMMANDER", "GENERAL"];
		const forward = shuffle([...free].filter((index) => !rear.has(rowOf(index))));
		if (forward.length < leaders.length) return null;
		leaders.forEach((kind, offset) => put(forward[offset], kind));
		const placed = board.reduce((counts, piece) => {
			if (piece?.owner === owner) counts[piece.kind] += 1;
			return counts;
		}, Object.fromEntries(KINDS.map((kind) => [kind, 0])));
		const remaining = shuffle(KINDS.flatMap((kind) => Array(COUNTS[kind] - placed[kind]).fill(kind)));
		const leftover = [...free];
		if (leftover.length !== remaining.length) return null;
		shuffle(leftover).forEach((index, offset) => put(index, remaining[offset]));
		return board;
	}
	function validateDeployment(board, owner) {
		const pieces = board.flatMap((piece, index) => piece?.owner === owner ? [{
			piece,
			index
		}] : []);
		if (pieces.length !== 25) return false;
		return pieces.every(({ piece, index }) => {
			if (piece.kind === "FLAG") return ownHeadquarters(owner).includes(index);
			if (piece.kind === "MINE") return rearRows(owner).has(rowOf(index));
			if (piece.kind === "BOMB") return rowOf(index) !== frontRow(owner);
			return deploymentSquares(owner).includes(index);
		});
	}
	//#endregion
	//#region lib/game.ts
	function newGame() {
		const bot = strategicDeployment("bot");
		const human = randomDeployment("human");
		return {
			board: bot.map((piece, index) => piece ?? human[index]),
			turn: "human",
			phase: "deployment",
			winner: null,
			moveCount: 0,
			history: [],
			knowledge: {},
			lastBotHint: null
		};
	}
	function startGame(state) {
		if (!validateDeployment(state.board, "human")) return state;
		return {
			...state,
			phase: "playing",
			turn: "human",
			opening: state.board.map((piece) => piece ? { ...piece } : null)
		};
	}
	function battleOutcome(attacker, defender) {
		if (defender === "FLAG") return 1;
		if (attacker === "BOMB" || defender === "BOMB") return 0;
		if (defender === "MINE") return attacker === "ENGINEER" ? 1 : -1;
		if (attacker === defender) return 0;
		return RANK[attacker] < RANK[defender] ? 1 : -1;
	}
	function revealFlag(board, owner) {
		const position = ownHeadquarters(owner).find((index) => board[index]?.owner === owner && board[index]?.kind === "FLAG");
		if (position === void 0) return null;
		board[position] = {
			...board[position],
			revealed: true
		};
		return position;
	}
	function candidatesThatBeat(kind, attacking) {
		return KINDS.filter((candidate) => {
			if (candidate === "FLAG" || candidate === "MINE") return false;
			const result = attacking ? battleOutcome(candidate, kind) : battleOutcome(kind, candidate);
			return result > 0 === attacking && result !== 0;
		});
	}
	function updateKnowledge(state, move, attacker, defender, outcome) {
		const knowledge = { ...state.knowledge };
		if (attacker.owner === "human") {
			const previous = knowledge[move.from] ?? KINDS.filter((kind) => kind !== "FLAG" && kind !== "MINE");
			delete knowledge[move.from];
			if (!defender) knowledge[move.to] = previous;
			else if (outcome === 1 && defender.owner === "bot") knowledge[move.to] = previous.filter((kind) => battleOutcome(kind, defender.kind) > 0);
			else delete knowledge[move.to];
		} else if (defender?.owner === "human") if (outcome === -1) {
			const previous = knowledge[move.to] ?? KINDS;
			knowledge[move.to] = previous.filter((kind) => candidatesThatBeat(attacker.kind, false).includes(kind));
		} else delete knowledge[move.to];
		return knowledge;
	}
	function applyMove(state, move, learningKey) {
		if (state.phase !== "playing") return state;
		if (!legalMoves(state.board, state.turn).some((candidate) => candidate.from === move.from && candidate.to === move.to)) return state;
		const board = [...state.board];
		const attacker = board[move.from];
		const defender = board[move.to];
		board[move.from] = null;
		let outcome = null;
		let winner = null;
		let message = `${coordinate(move.from)} → ${coordinate(move.to)}`;
		const revealed = [];
		if (!defender) board[move.to] = attacker;
		else {
			outcome = battleOutcome(attacker.kind, defender.kind);
			if (attacker.kind === "COMMANDER" && outcome <= 0) {
				const flag = revealFlag(board, attacker.owner);
				if (flag !== null) revealed.push(flag);
			}
			if (defender.kind === "COMMANDER" && outcome >= 0) {
				const flag = revealFlag(board, defender.owner);
				if (flag !== null) revealed.push(flag);
			}
			board[move.to] = outcome > 0 ? {
				...attacker,
				revealed: false
			} : outcome < 0 ? {
				...defender,
				revealed: false
			} : null;
			if (defender.kind === "FLAG") {
				winner = attacker.owner;
				message += " · 军旗被夺";
			} else message += outcome > 0 ? " · 攻方胜" : outcome < 0 ? " · 守方胜" : " · 同归于尽";
			if (revealed.length) message += ` · 亮旗 ${revealed.map(coordinate).join(" / ")}`;
		}
		const knowledge = updateKnowledge(state, move, attacker, defender, outcome);
		const next = other(state.turn);
		if (!winner && legalMoves(board, next).length === 0) winner = state.turn;
		const humanPieceDied = attacker.owner === "bot" && defender?.owner === "human" && outcome !== null && outcome >= 0;
		const lastBotHint = attacker.owner === "bot" ? {
			square: move.to,
			victimKind: humanPieceDied ? defender.kind : void 0,
			sequence: state.moveCount + 1
		} : state.lastBotHint;
		return {
			...state,
			board,
			turn: next,
			phase: winner ? "finished" : "playing",
			winner,
			moveCount: state.moveCount + 1,
			history: [...state.history, {
				move,
				owner: attacker.owner,
				message,
				learningKey
			}],
			knowledge,
			lastBotHint
		};
	}
	//#endregion
	//#region lib/memory.ts
	function emptyMemory() {
		return {
			games: 0,
			actionValues: {},
			matches: []
		};
	}
	function learningKey(state, move) {
		const actor = state.board[move.from];
		const target = state.board[move.to];
		const targetClass = !target ? "empty" : target.revealed && target.kind === "FLAG" ? "flag" : "unknown";
		const phaseBucket = Math.min(4, Math.floor(state.moveCount / 12));
		return `${actor?.kind ?? "?"}|${move.from}-${move.to}|${targetClass}|${phaseBucket}`;
	}
	function learnedBonus(memory, key) {
		const value = memory.actionValues[key];
		if (!value || value.visits < 1) return 0;
		const confidence = Math.min(1, value.visits / 4);
		return value.total / value.visits * 7 * confidence;
	}
	//#endregion
	//#region lib/weights.ts
	var WEIGHTS = {
		capture: 2.8,
		flag_capture: 137.5478244461,
		forward: .4070295049,
		camp: 1.1,
		mobility: .0490267001,
		protect_flag: .4024492383,
		revealed_flag_hunt: 4.2620660189,
		unknown_risk: .0907965881,
		belief_battle: 9.5539998758,
		hq_pressure: .9435573435,
		hq_strike: 15.7040728292,
		mine_risk: .5,
		engineer_mine: 3.9957507745,
		engineer_waste: -12,
		noise: .1594516286,
		eval_material: 1.7,
		eval_mobility: .0522149392,
		eval_terminal: 2e3,
		eval_hq_attack: .9,
		eval_hq_attack_certain: 3.0905818961,
		eval_hq_defense: .8596938267,
		eval_hq_defense_certain: 3.0566762473,
		eval_hq_breach: 26,
		eval_hq_guard: 5.5,
		eval_commander: 0,
		engineer_expose: -6,
		eval_immobilize: 110
	};
	//#endregion
	//#region lib/bot.ts
	var EVAL_HORIZON = 6;
	var SETTINGS = {
		casual: {
			beam: 5,
			continuations: 0,
			samples: 4,
			noise: 3.5,
			replies: 1
		},
		focused: {
			beam: 12,
			continuations: 1,
			samples: 14,
			noise: .7,
			replies: 4
		},
		deep: {
			beam: 18,
			continuations: 2,
			samples: 28,
			noise: .12,
			replies: 5
		}
	};
	function distance(left, right) {
		return Math.abs(rowOf(left) - rowOf(right)) + Math.abs(colOf(left) - colOf(right));
	}
	function pieceValue(kind) {
		if (kind === "FLAG") return 50;
		if (kind === "BOMB") return 7;
		if (kind === "MINE") return 5;
		return 12 - RANK[kind];
	}
	function expectedBattle(attacker, possible) {
		if (!possible.length) return 0;
		return possible.reduce((sum, kind) => sum + battleOutcome(attacker, kind), 0) / possible.length;
	}
	function nearest(from, targets) {
		return targets.reduce((best, target) => Math.min(best, distance(from, target)), Number.POSITIVE_INFINITY);
	}
	function scoreMove(state, move, owner) {
		const piece = state.board[move.from];
		const target = state.board[move.to];
		let score = (owner === "bot" ? 1 : -1) * (rowOf(move.to) - rowOf(move.from)) * WEIGHTS.forward;
		if (CAMPS.has(move.to)) score += WEIGHTS.camp;
		if (piece.kind === "ENGINEER" && revealsEngineer(state.board, move.from, move.to)) score += WEIGHTS.engineer_expose;
		const enemy = other(owner);
		const enemyFlagSquares = liveFlagSquares(state.board, enemy, enemy === "human" ? state.knowledge : void 0);
		const certain = enemyFlagSquares.length === 1;
		if (target) {
			score += WEIGHTS.capture;
			if (target.revealed || enemyFlagSquares.includes(move.to)) score += target.revealed || certain ? WEIGHTS.flag_capture : WEIGHTS.hq_strike;
			else if (owner === "bot" && state.knowledge[move.to]) score += expectedBattle(piece.kind, state.knowledge[move.to]) * WEIGHTS.belief_battle;
			else if (rearRows(enemy).has(rowOf(move.to))) score += piece.kind === "ENGINEER" ? WEIGHTS.engineer_mine : -pieceValue(piece.kind) * WEIGHTS.mine_risk;
			else if (piece.kind === "ENGINEER") score += WEIGHTS.engineer_waste;
			else score -= pieceValue(piece.kind) * WEIGHTS.unknown_risk;
		}
		const ownFlagSquares = liveFlagSquares(state.board, owner);
		if (ownFlagSquares.length) score += (nearest(move.from, ownFlagSquares) - nearest(move.to, ownFlagSquares)) * WEIGHTS.protect_flag;
		if (enemyFlagSquares.length) score += (nearest(move.from, enemyFlagSquares) - nearest(move.to, enemyFlagSquares)) * (certain ? WEIGHTS.revealed_flag_hunt : WEIGHTS.hq_pressure);
		return score;
	}
	/** Distance in *moves*, so a raider sitting on a railway counts as near. */
	function closestRaider(board, side, targets) {
		if (!targets.length) return null;
		let best = null;
		board.forEach((piece, index) => {
			if (!piece || piece.owner !== side) return;
			if (piece.kind === "FLAG" || piece.kind === "MINE") return;
			if (HEADQUARTERS.has(index)) return;
			const reach = targets.reduce((closest, target) => Math.min(closest, moveDistance(index, target)), Number.POSITIVE_INFINITY);
			if (best === null || reach < best) best = reach;
		});
		return best;
	}
	/**
	* Reward closing on the enemy headquarters and punish losing our own. Without
	* this the evaluation is pure material and the bot never plays for the win.
	*/
	function flagPressure(state) {
		let value = 0;
		const attacking = liveFlagSquares(state.board, "human", state.knowledge);
		const reach = closestRaider(state.board, "bot", attacking);
		if (reach !== null) value += Math.max(0, EVAL_HORIZON - reach) * (attacking.length === 1 ? WEIGHTS.eval_hq_attack_certain : WEIGHTS.eval_hq_attack);
		const defending = liveFlagSquares(state.board, "bot");
		const threat = closestRaider(state.board, "human", defending);
		if (threat !== null) {
			value -= Math.max(0, EVAL_HORIZON - threat) * (defending.length === 1 ? WEIGHTS.eval_hq_defense_certain : WEIGHTS.eval_hq_defense);
			if (threat <= 1) value -= WEIGHTS.eval_hq_breach;
		}
		value += WEIGHTS.eval_hq_guard * guards(state.board, "bot", defending);
		return value;
	}
	/** Own movable pieces shielding a headquarters that may hold our flag. */
	function guards(board, owner, squares) {
		if (!squares.length) return 0;
		let total = 0;
		board.forEach((piece, index) => {
			if (!piece || piece.owner !== owner) return;
			if (piece.kind === "FLAG" || piece.kind === "MINE") return;
			if (HEADQUARTERS.has(index)) return;
			if (squares.some((square) => distance(index, square) === 1)) total += 1;
		});
		return total;
	}
	function evaluateState(state) {
		if (state.winner === "bot") return WEIGHTS.eval_terminal;
		if (state.winner === "human") return -WEIGHTS.eval_terminal;
		const material = state.board.reduce((sum, piece) => {
			if (!piece) return sum;
			return sum + pieceValue(piece.kind) * (piece.owner === "bot" ? 1 : -1);
		}, 0);
		const mobility = legalMoves(state.board, "bot").length - legalMoves(state.board, "human").length;
		const concealment = commanderShield(state.board, "bot") - commanderShield(state.board, "human");
		const squeeze = WEIGHTS.eval_immobilize * (2 ** -mobileCount(state.board, "human") - 2 ** -mobileCount(state.board, "bot"));
		return material * WEIGHTS.eval_material + mobility * WEIGHTS.eval_mobility + concealment * WEIGHTS.eval_commander + squeeze + flagPressure(state);
	}
	/** Pieces this side can actually move. */
	function mobileCount(board, owner) {
		let total = 0;
		board.forEach((piece, index) => {
			if (!piece || piece.owner !== owner) return;
			if (piece.kind === "FLAG" || piece.kind === "MINE") return;
			if (HEADQUARTERS.has(index)) return;
			total += 1;
		});
		return total;
	}
	/**
	* 1 while this side's commander lives and its flag is still hidden.
	*
	* The commander's life is worth exactly the concealment it buys, so the term
	* vanishes once the flag is out. Evaluated on sampled worlds, so reading the
	* opponent's commander is a guess, not a peek.
	*/
	function commanderShield(board, owner) {
		const squares = liveFlagSquares(board, owner);
		if (!squares.length || squares.some((index) => board[index]?.revealed)) return 0;
		return board.some((piece) => piece?.owner === owner && piece.kind === "COMMANDER") ? 1 : 0;
	}
	function continuationValue(state, width) {
		if (state.phase === "finished" || width === 0) return evaluateState(state);
		const continuations = legalMoves(state.board, "bot").map((move) => ({
			move,
			score: scoreMove(state, move, "bot")
		})).sort((left, right) => right.score - left.score).slice(0, width);
		if (!continuations.length) return evaluateState(state);
		return Math.max(...continuations.map(({ move }) => evaluateState(applyMove(state, move))));
	}
	function searchValue(state, move, replyWidth, continuationWidth) {
		const next = applyMove(state, move);
		if (next.phase === "finished") return evaluateState(next);
		const replies = legalMoves(next.board, "human").map((reply) => ({
			reply,
			score: scoreMove(next, reply, "human")
		})).sort((left, right) => right.score - left.score).slice(0, replyWidth);
		if (!replies.length) return evaluateState(next);
		return Math.min(...replies.map(({ reply }) => continuationValue(applyMove(next, reply), continuationWidth)));
	}
	function chooseBotMove(state, difficulty, memory) {
		const settings = SETTINGS[difficulty];
		const ranked = legalMoves(state.board, "bot").map((move) => ({
			move,
			base: scoreMove(state, move, "bot") + learnedBonus(memory, learningKey(state, move))
		})).sort((left, right) => right.base - left.base).slice(0, settings.beam);
		if (!ranked.length) throw new Error("The bot has no legal move to choose.");
		const worlds = sampleWorlds(state, settings.samples);
		return ranked.map(({ move, base }) => {
			const search = worlds.length ? worlds.reduce((sum, world) => sum + searchValue(world, move, settings.replies, settings.continuations), 0) / worlds.length : 0;
			return {
				move,
				score: base * 2 + search + (Math.random() - .5) * settings.noise
			};
		}).reduce((best, candidate) => candidate.score > best.score ? candidate : best).move;
	}
	//#endregion
	//#region lib/engine-api.ts
	/**
	* A standalone facade for driving the engine from outside this app.
	*
	* Built to `dist-engine/junqi-engine.js` by `npm run build:engine`, which can
	* be pasted into any page's console. The point is that the adapter for a
	* foreign board only has to describe what it can *see* -- this engine is an
	* information-set player and never reads a hidden rank, so "occupancy plus my
	* own pieces" is exactly its natural input.
	*
	* Coordinates are the usual A1..E12, A1 being the north-west corner (the bot's
	* side). If the foreign board is oriented the other way, flip it in the
	* adapter, not here.
	*/
	function parseCoordinate(square) {
		const file = square.trim().toUpperCase();
		const col = file.charCodeAt(0) - 65;
		const row = Number(file.slice(1)) - 1;
		if (!(col >= 0 && col < 5) || !(row >= 0 && row < 12)) throw new Error(`not a square on a 12x5 board: ${square}`);
		return indexOf(row, col);
	}
	/**
	* An enemy piece whose rank we cannot see still needs *a* rank internally.
	* Anything unrevealed is overwritten by the belief sampler before it is
	* searched, so the placeholder never reaches the evaluation.
	*/
	var PLACEHOLDER = "CAPTAIN";
	function toGameState(external, turn = "bot") {
		const board = Array(60).fill(null);
		for (const [square, cell] of Object.entries(external)) {
			const index = parseCoordinate(square);
			board[index] = {
				owner: cell.owner,
				kind: cell.kind ?? PLACEHOLDER,
				revealed: Boolean(cell.revealed)
			};
		}
		return {
			board,
			turn,
			phase: "playing",
			winner: null,
			moveCount: 0,
			history: [],
			knowledge: {},
			lastBotHint: null,
			opening: board.map((piece) => piece ? { ...piece } : null)
		};
	}
	/**
	* Ranks that survive an attack by `attacker`.
	*
	* What you learn when your own piece dies attacking a square. Note FLAG is
	* absent -- a flag loses to everything -- so a failed attack on a headquarters
	* proves that square is not the flag.
	*/
	function survivorsOf(attacker) {
		return KINDS.filter((defender) => battleOutcome(attacker, defender) < 0);
	}
	/**
	* Pick a move for "bot" from a board described in A1..E12 coordinates.
	*
	* `knowledge` carries what past battles proved about enemy squares, keyed the
	* same way. Without it the engine re-attacks a square that already killed a
	* stronger piece of ours, because it has no way to remember that it did.
	*/
	function chooseMove(external, difficulty = "deep", knowledge = {}) {
		const state = toGameState(external, "bot");
		for (const [square, kinds] of Object.entries(knowledge)) state.knowledge[parseCoordinate(square)] = kinds;
		const move = chooseBotMove(state, difficulty, emptyMemory());
		return {
			from: coordinate(move.from),
			to: coordinate(move.to)
		};
	}
	/** A fresh legal opening for our side, for adapters that must deploy first. */
	function suggestDeployment() {
		const state = startGame(newGame());
		const layout = {};
		state.board.forEach((piece, index) => {
			if (piece?.owner === "bot") layout[coordinate(index)] = piece.kind;
		});
		return layout;
	}
	//#endregion
	exports.COLS = COLS;
	exports.ROWS = ROWS;
	exports.chooseMove = chooseMove;
	exports.colOf = colOf;
	exports.coordinate = coordinate;
	exports.indexOf = indexOf;
	exports.rowOf = rowOf;
	exports.suggestDeployment = suggestDeployment;
	exports.survivorsOf = survivorsOf;
	exports.toGameState = toGameState;
	return exports;
})({});
