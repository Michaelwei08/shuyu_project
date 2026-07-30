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


def road_neighbors(position: Position) -> Iterable[Position]:
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


def railway_neighbors(position: Position) -> Iterable[Position]:
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
    if source not in RAILWAYS:
        return set()
    row, column = source
    destinations: set[Position] = set()
    for row_step, column_step in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        cursor = (row + row_step, column + column_step)
        while cursor in RAILWAYS:
            crosses_river = {cursor[0], cursor[0] - row_step} == {5, 6}
            if crosses_river and cursor[1] not in (0, 2, 4):
                break
            destinations.add(cursor)
            if cursor in occupied:
                break
            cursor = (cursor[0] + row_step, cursor[1] + column_step)
    return destinations
