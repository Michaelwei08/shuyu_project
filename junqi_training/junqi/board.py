from __future__ import annotations

from collections import deque
from collections.abc import Iterable

from .types import Position

ROWS = 12
COLUMNS = 5

CAMPS: frozenset[Position] = frozenset(
    {
        (2, 1),
        (2, 3),
        (3, 2),
        (4, 1),
        (4, 3),
        (7, 1),
        (7, 3),
        (8, 2),
        (9, 1),
        (9, 3),
    }
)
HEADQUARTERS: frozenset[Position] = frozenset(
    {(0, 1), (0, 3), (11, 1), (11, 3)}
)
RAILWAYS: frozenset[Position] = frozenset(
    {(row, column) for row in (1, 5, 6, 10) for column in range(COLUMNS)}
    | {(row, column) for row in range(1, 11) for column in (0, 4)}
    | {(5, 2), (6, 2)}
)


def in_bounds(position: Position) -> bool:
    row, column = position
    return 0 <= row < ROWS and 0 <= column < COLUMNS


def _road_neighbors(position: Position) -> Iterable[Position]:
    row, column = position
    offsets = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    if position in CAMPS:
        offsets += [(1, 1), (1, -1), (-1, 1), (-1, -1)]
    for offset_row, offset_column in offsets:
        candidate = (row + offset_row, column + offset_column)
        crosses_river = {row, candidate[0]} == {5, 6}
        if in_bounds(candidate) and not (crosses_river and column not in (0, 2, 4)):
            yield candidate

    for camp in CAMPS:
        if abs(row - camp[0]) == 1 and abs(column - camp[1]) == 1:
            yield camp


def _railway_neighbors(position: Position) -> Iterable[Position]:
    if position not in RAILWAYS:
        return
    row, column = position
    for candidate in (
        (row + 1, column),
        (row - 1, column),
        (row, column + 1),
        (row, column - 1),
    ):
        crosses_river = {row, candidate[0]} == {5, 6}
        if candidate in RAILWAYS and not (
            crosses_river and column not in (0, 2, 4)
        ):
            yield candidate


def _rail_rays(position: Position) -> tuple[tuple[Position, ...], ...]:
    """Squares reachable in a straight rail slide, before blocking is applied."""
    if position not in RAILWAYS:
        return ()
    row, column = position
    rays: list[tuple[Position, ...]] = []
    for row_step, column_step in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        ray: list[Position] = []
        cursor = (row + row_step, column + column_step)
        while cursor in RAILWAYS:
            crosses_river = {cursor[0], cursor[0] - row_step} == {5, 6}
            if crosses_river and cursor[1] not in (0, 2, 4):
                break
            ray.append(cursor)
            cursor = (cursor[0] + row_step, cursor[1] + column_step)
        if ray:
            rays.append(tuple(ray))
    return tuple(rays)


# The board graph is static, so build it once. Recomputing `road_neighbors` per
# call cost roughly a third of all search time -- it was being called about
# three million times per dozen plies.
_SQUARES = [(row, column) for row in range(ROWS) for column in range(COLUMNS)]
ROAD_NEIGHBORS: dict[Position, frozenset[Position]] = {
    square: frozenset(_road_neighbors(square)) for square in _SQUARES
}
RAILWAY_NEIGHBORS: dict[Position, tuple[Position, ...]] = {
    square: tuple(_railway_neighbors(square)) for square in _SQUARES
}
RAIL_RAYS: dict[Position, tuple[tuple[Position, ...], ...]] = {
    square: _rail_rays(square) for square in _SQUARES
}


def road_neighbors(position: Position) -> frozenset[Position]:
    return ROAD_NEIGHBORS[position]


def railway_neighbors(position: Position) -> tuple[Position, ...]:
    return RAILWAY_NEIGHBORS[position]


def engineer_rail_destinations(
    source: Position, occupied: set[Position]
) -> set[Position]:
    """Return reachable rail squares; occupied squares are endpoints, not transit."""
    if source not in RAILWAYS:
        return set()
    destinations: set[Position] = set()
    visited = {source}
    queue = deque([source])
    while queue:
        current = queue.popleft()
        for neighbor in railway_neighbors(current):
            if neighbor in visited:
                continue
            visited.add(neighbor)
            destinations.add(neighbor)
            if neighbor not in occupied:
                queue.append(neighbor)
    return destinations


def straight_rail_destinations(
    source: Position, occupied: set[Position]
) -> set[Position]:
    destinations: set[Position] = set()
    for ray in RAIL_RAYS[source]:
        for square in ray:
            destinations.add(square)
            if square in occupied:
                break
    return destinations


def _move_distances() -> dict[Position, dict[Position, int]]:
    """Moves -- not squares -- between every pair, on an empty board.

    Manhattan distance is badly wrong here: a piece on a railway crosses the
    whole rank in one move, so E2 is two moves from B1 while being Manhattan 4
    away. Threat and pressure both have to be counted in moves or the bot never
    sees a rail raid coming.

    Blockers only ever slow a piece down, so the empty board gives a lower
    bound -- the conservative direction for judging danger.
    """
    table: dict[Position, dict[Position, int]] = {}
    empty: set[Position] = set()
    for source in _SQUARES:
        distances = {source: 0}
        queue = deque([source])
        while queue:
            current = queue.popleft()
            step = distances[current] + 1
            for neighbor in (
                ROAD_NEIGHBORS[current] | straight_rail_destinations(current, empty)
            ):
                if neighbor not in distances:
                    distances[neighbor] = step
                    queue.append(neighbor)
        table[source] = distances
    return table


MOVE_DISTANCE: dict[Position, dict[Position, int]] = _move_distances()
UNREACHABLE = ROWS * COLUMNS


def move_distance(source: Position, target: Position) -> int:
    return MOVE_DISTANCE[source].get(target, UNREACHABLE)
