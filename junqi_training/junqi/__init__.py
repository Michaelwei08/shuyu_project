"""Trainable Junqi game prototype."""

import sys

if sys.version_info < (3, 10):
    # Fail with the actual problem rather than a TypeError about `slots` from
    # deep inside types.py. This bites on shared machines whose system python
    # is old and whose conda env is easy to lose between shells.
    raise RuntimeError(
        "junqi needs Python 3.10 or newer (dataclass slots, zip strict), but "
        f"this is {sys.version.split()[0]} at {sys.executable}.\n"
        "Activate the right environment, or call its interpreter directly:\n"
        "    conda activate junqi\n"
        "    ~/.conda/envs/junqi/bin/python -m junqi.<module>"
    )

from .game import Game  # noqa: E402  - must follow the version guard
from .types import Move, Owner, Piece, PieceKind  # noqa: E402

__all__ = ["Game", "Move", "Owner", "Piece", "PieceKind"]
