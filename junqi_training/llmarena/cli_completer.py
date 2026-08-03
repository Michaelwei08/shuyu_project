"""Drive the model through the `claude` CLI instead of the Messages API.

Runs on whatever credentials Claude Code already has -- including a
subscription -- so it costs no API spend. That makes it the right instrument
for a cheap first look, and the wrong one for the headline number. The reason
is measured, not asserted:

    plain `claude -p`             cache_creation_input_tokens = 25,830
    with --strict-mcp-config      cache_creation_input_tokens =  5,153

Even at its leanest there are ~5.1k tokens of Claude Code harness -- its own
system prompt, its built-in tool definitions -- between the experiment and the
model, against a ~1.6k-token experiment prompt. It is constant across
conditions, so a *difference* between two scaffolds is still roughly
interpretable; it is not constant across harnesses, so the absolute
search-equivalent number is not comparable to a clean API measurement or to
anyone else's. Use this to decide whether the result is interesting. Pay for
the number you intend to publish.

Two hazards, both guarded below:

* **CLAUDE.md auto-discovery.** Non-`--bare` Claude Code loads the CLAUDE.md of
  its working directory. Run this from `military/` and the model is handed the
  project's own strategy documentation -- the mine prior, the flag-candidate
  deduction, the blind-attack pricing. That is not a leak of hidden ranks, but
  it is a total giveaway of the answers the probes are asking for.
  :func:`claude_cli_completer` refuses to start in a directory containing one.
* **`--bare` does not help.** It would strip the harness, but its auth is
  strictly `ANTHROPIC_API_KEY` -- OAuth and keychain are never read -- so it
  defeats the entire point of using the subscription.
"""

from __future__ import annotations

import json
import subprocess
from collections.abc import Callable
from pathlib import Path

from .anthropic_completer import Usage, split_cacheable


class HarnessContamination(RuntimeError):
    """Raised when the working directory would inject project context."""


def _check_clean(workdir: Path) -> None:
    for name in ("CLAUDE.md", "AGENTS.md", ".claude"):
        if (workdir / name).exists():
            raise HarnessContamination(
                f"{workdir} contains {name}; Claude Code would load it into the "
                "prompt and hand the model the answers. Point --workdir at an "
                "empty scratch directory."
            )


def claude_cli_completer(
    workdir: str | Path,
    *,
    model: str | None = None,
    timeout: float = 180.0,
    retries: int = 2,
    usage: Usage | None = None,
    executable: str = "claude",
) -> Callable[[str], str]:
    """Build a ``(prompt) -> text`` completer backed by ``claude -p``.

    The rules block goes to ``--system-prompt`` and the position to stdin, which
    mirrors the API path's system/user split so the two are comparable to each
    other even though neither is comparable to a bare model.
    """
    directory = Path(workdir)
    directory.mkdir(parents=True, exist_ok=True)
    _check_clean(directory)

    # `--strict-mcp-config` on its own is not enough, and the failure is
    # expensive rather than loud. With the flag alone the harness prefix was
    # measured at 5,159 tokens on one call and 25,828 on the next -- i.e. the
    # session's MCP servers sometimes load anyway. When they do, the CLI also
    # tries to *connect* to each one, which is the most likely reason a whole
    # match run timed out at 180s per call while a bare prompt answered in six
    # seconds. Pointing --mcp-config at an empty server set makes it
    # deterministic: 5,159 every time.
    # Absolute: the subprocess runs with `cwd=directory`, so a relative path
    # would resolve against it a second time.
    mcp_config = (directory / "empty-mcp.json").resolve()
    mcp_config.write_text('{"mcpServers":{}}', encoding="utf-8")

    def complete(prompt: str) -> str:
        system, rest = split_cacheable(prompt)
        command = [
            executable,
            "--print",
            "--output-format",
            "json",
            # Both are needed. Measured: 25,830 -> 5,159, and deterministic.
            "--mcp-config",
            str(mcp_config),
            "--strict-mcp-config",
            "--no-session-persistence",
        ]
        if system is not None:
            command += ["--system-prompt", system]
        if model is not None:
            command += ["--model", model]

        # A stalled call kills the whole game, and `play_match` cannot resume
        # mid-game -- only a re-run can, off the cache. Measured latency is
        # ~17-23s against a 180s ceiling, so a timeout here is a transient
        # stall rather than genuine slowness, and retrying is far cheaper than
        # losing the game.
        for attempt in range(retries + 1):
            try:
                finished = subprocess.run(
                    command,
                    input=rest,
                    capture_output=True,
                    text=True,
                    encoding="utf-8",
                    cwd=directory,
                    timeout=timeout,
                )
                break
            except subprocess.TimeoutExpired:
                if attempt == retries:
                    raise
        if finished.returncode != 0:
            raise RuntimeError(
                f"claude exited {finished.returncode}: {finished.stderr[:400]}"
            )
        payload = json.loads(finished.stdout)
        if payload.get("is_error"):
            raise RuntimeError(f"claude reported an error: {payload.get('result')}")

        if usage is not None:
            counts = payload.get("usage", {})
            usage.record(
                _Counts(
                    input_tokens=counts.get("input_tokens", 0),
                    output_tokens=counts.get("output_tokens", 0),
                    cache_read_input_tokens=counts.get("cache_read_input_tokens", 0),
                    cache_creation_input_tokens=counts.get(
                        "cache_creation_input_tokens", 0
                    ),
                )
            )
        return payload.get("result", "")

    return complete


class _Counts:
    """Duck-types the SDK usage object that :meth:`Usage.record` reads."""

    __slots__ = (
        "input_tokens",
        "output_tokens",
        "cache_read_input_tokens",
        "cache_creation_input_tokens",
    )

    def __init__(self, **fields: int) -> None:
        for key, value in fields.items():
            setattr(self, key, value)


def register() -> None:
    """Add ``claude-cli`` to the completer registry used by ``AgentSpec``."""
    from .agent import COMPLETERS

    COMPLETERS["claude-cli"] = lambda options: claude_cli_completer(
        options.get("workdir", "data/cli-scratch"),
        model=options.get("model") or None,
        timeout=float(options.get("timeout", "180")),
        retries=int(options.get("retries", "2")),
    )
