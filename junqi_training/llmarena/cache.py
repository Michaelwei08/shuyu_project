"""A persistent prompt -> response cache.

One cheap component that buys four separate things, which is why it goes in on
day one rather than when the bill arrives:

* **Reproducibility.** The rest of the harness is built on common random
  numbers -- ``arena.compare`` only means anything because two candidates meet
  byte-identical states. An LLM breaks that even at temperature 0, since
  providers do not guarantee identical samples. Replaying through a warm cache
  restores exact reproducibility for everything downstream of the call.
* **Crash resume.** A game is ~42 calls and a sweep is thousands; a run that
  dies at 80% must not re-spend the first 80%.
* **Offline re-scoring.** Raw responses are kept verbatim, so a change to the
  move parser or the repair policy is re-scored without re-spending anything.
* **Cost.** Repeated positions -- openings above all -- hit.

Layout is one file per entry, sharded two hex characters deep, written to a
temporary name and then ``os.replace``d. That is atomic on both NTFS and POSIX,
so the process pool can write concurrently without a lock. A single append-only
JSONL would have been smaller and would have interleaved corrupt lines under
``ProcessPoolExecutor``.
"""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


def cache_key(model: str, scaffold: str, prompt: str, variant: int = 0) -> str:
    """Hash the full call identity.

    ``variant`` separates deliberate re-asks of an identical prompt -- a retry
    after an illegal move must be allowed to draw a fresh sample instead of
    replaying the rejected one out of the cache.
    """
    digest = hashlib.sha256()
    for part in (model, scaffold, str(variant), prompt):
        digest.update(part.encode("utf-8"))
        digest.update(b"\x00")
    return digest.hexdigest()


@dataclass
class PromptCache:
    """Content-addressed store of completions."""

    root: Path
    hits: int = field(default=0, init=False)
    misses: int = field(default=0, init=False)
    writes: int = field(default=0, init=False)

    def __post_init__(self) -> None:
        self.root = Path(self.root)

    def _path(self, key: str) -> Path:
        return self.root / key[:2] / f"{key}.json"

    def get(
        self, model: str, scaffold: str, prompt: str, variant: int = 0
    ) -> str | None:
        path = self._path(cache_key(model, scaffold, prompt, variant))
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (FileNotFoundError, json.JSONDecodeError):
            self.misses += 1
            return None
        self.hits += 1
        return payload["response"]

    def put(
        self,
        model: str,
        scaffold: str,
        prompt: str,
        response: str,
        variant: int = 0,
        meta: dict[str, Any] | None = None,
    ) -> None:
        key = cache_key(model, scaffold, prompt, variant)
        path = self._path(key)
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "model": model,
            "scaffold": scaffold,
            "variant": variant,
            "prompt": prompt,
            "response": response,
            "meta": meta or {},
        }
        # Write-then-rename: a reader never sees a half-written entry, and two
        # workers racing on the same key both end up with a valid file.
        handle, temporary = tempfile.mkstemp(dir=path.parent, suffix=".tmp")
        try:
            with os.fdopen(handle, "w", encoding="utf-8") as stream:
                json.dump(payload, stream, ensure_ascii=False)
            os.replace(temporary, path)
        except BaseException:
            Path(temporary).unlink(missing_ok=True)
            raise
        self.writes += 1

    def stats(self) -> dict[str, int]:
        return {"hits": self.hits, "misses": self.misses, "writes": self.writes}


class NullCache:
    """Drop-in that stores nothing, for tests that want every call to land."""

    def get(self, *_args: Any, **_kwargs: Any) -> str | None:
        return None

    def put(self, *_args: Any, **_kwargs: Any) -> None:
        return None

    def stats(self) -> dict[str, int]:
        return {"hits": 0, "misses": 0, "writes": 0}
