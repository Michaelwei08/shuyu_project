"""LLM agents against the Junqi engine.

Kept out of ``junqi`` on purpose: that package is stdlib-only and the rules
engine should stay that way. Nothing here is imported by the engine; the arena
reaches an agent in this package through
``AgentSpec(builder="llmarena.agent:build_agent")``, which resolves by import
inside the worker process.

The load-bearing piece is :mod:`llmarena.view` -- see its docstring before
touching anything that produces prompt text.
"""

from .agent import LLMAgent, TurnRecord, first_illegal_ply, parse_response
from .belief import BeliefTracker
from .cache import NullCache, PromptCache
from .view import SCAFFOLDS, Observation, Scaffold, build_observation, render

__all__ = [
    "BeliefTracker",
    "LLMAgent",
    "NullCache",
    "Observation",
    "PromptCache",
    "SCAFFOLDS",
    "Scaffold",
    "TurnRecord",
    "build_observation",
    "first_illegal_ply",
    "parse_response",
    "render",
]
