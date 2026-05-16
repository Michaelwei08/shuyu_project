from __future__ import annotations

import json
import re
from collections import Counter
from html import escape
from pathlib import Path
from typing import Any
import zipfile
from xml.etree import ElementTree as ET

import yaml


PROJECT_DIR = Path(__file__).resolve().parent
OUTLINE_PATH = PROJECT_DIR / "outline.yaml"
FIELDS_PATH = PROJECT_DIR / "fields.yaml"
OUTPUT_PATH = PROJECT_DIR / "report.md"
VISUALS_DIR = PROJECT_DIR / "visuals"
DATA_DIR = PROJECT_DIR / "data"
SOURCE_SUMMARY_PATH = PROJECT_DIR / "validation" / "source_summary.yaml"
SOURCE_VALIDATION_JSON_PATH = PROJECT_DIR / "validation" / "source_validation.json"
CLAIM_VALIDATION_PATH = PROJECT_DIR / "validation" / "claim_validation.md"
BENCHMARK_DOCX_PATH = (
    PROJECT_DIR
    / "external_data"
    / "tools_assessment_hervk_SR-WGS"
    / "supplementary_table_3.docx"
)
BENCHMARK_CSV_PATH = DATA_DIR / "herv_tool_benchmark_table_3.csv"

CATEGORY_MAPPING = {
    "Basic Info": ["basic_info", "Basic Info"],
    "Technical Features": ["technical_features", "technical_characteristics", "Technical Features"],
    "Performance Metrics": ["performance_metrics", "performance", "Performance Metrics"],
    "Milestone Significance": ["milestone_significance", "milestones", "Milestone Significance"],
    "Business Info": ["business_info", "commercial_info", "Business Info"],
    "Competition & Ecosystem": ["competition_ecosystem", "competition", "Competition & Ecosystem"],
    "History": ["history", "History"],
    "Market Positioning": ["market_positioning", "market", "Market Positioning"],
}

INTERNAL_KEYS = {"_source_file", "uncertain"}
SUMMARY_FIELDS = ["category", "recommended_role"]


def load_yaml(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def slugify(value: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", value.lower()).strip("-")
    return slug or "item"


def local_md_link(path: Path) -> str:
    return f"[{path.name}](/{path.as_posix()})"


def local_image(path: Path, alt: str) -> str:
    relative = path.relative_to(PROJECT_DIR).as_posix()
    return f"![{alt}]({relative})"


def category_aliases() -> set[str]:
    aliases: set[str] = set()
    for values in CATEGORY_MAPPING.values():
        aliases.update(values)
    return aliases


def collect_all_fields(data: Any, field_names: set[str]) -> dict[str, Any]:
    found: dict[str, Any] = {}
    aliases = category_aliases()

    def visit(node: Any) -> None:
        if isinstance(node, dict):
            for key, value in node.items():
                if key in INTERNAL_KEYS:
                    continue
                if key in field_names:
                    found.setdefault(key, value)
                elif key in aliases and isinstance(value, dict):
                    visit(value)
                elif isinstance(value, (dict, list)):
                    visit(value)
        elif isinstance(node, list):
            for item in node:
                visit(item)

    visit(data)
    return found


def is_uncertain_value(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, str):
        return not value.strip() or "[uncertain]" in value
    if isinstance(value, list):
        if not value:
            return True
        return all(is_uncertain_value(item) for item in value)
    if isinstance(value, dict):
        if not value:
            return True
        return all(is_uncertain_value(item) for item in value.values())
    return False


def format_dict_inline(data: dict[str, Any]) -> str:
    parts = []
    for key, value in data.items():
        formatted = format_value(value, compact=True)
        if formatted:
            parts.append(f"{key}: {formatted}")
    return " | ".join(parts)


def format_list(values: list[Any], compact: bool = False) -> str:
    if not values:
        return ""
    if all(isinstance(item, dict) for item in values):
        lines = [format_dict_inline(item) for item in values if format_dict_inline(item)]
        if compact:
            return "; ".join(lines)
        return "\n".join(f"- {line}" for line in lines)
    if all(not isinstance(item, (dict, list)) for item in values):
        items = [str(item).strip() for item in values if str(item).strip()]
        if compact or sum(len(item) for item in items) <= 100:
            return ", ".join(items)
        return "\n".join(f"- {item}" for item in items)
    lines = [format_value(item, compact=compact) for item in values]
    lines = [line for line in lines if line]
    if compact:
        return "; ".join(lines)
    return "\n".join(f"- {line}" for line in lines)


def format_value(value: Any, compact: bool = False) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        text = value.strip()
        if not text or "[uncertain]" in text:
            return ""
        if compact:
            return re.sub(r"\s+", " ", text)
        if len(text) > 100:
            return "> " + text
        return text
    if isinstance(value, list):
        return format_list(value, compact=compact)
    if isinstance(value, dict):
        if compact:
            return format_dict_inline(value)
        lines = []
        for key, nested in value.items():
            formatted = format_value(nested, compact=True)
            if formatted:
                lines.append(f"- {key}: {formatted}")
        return "\n".join(lines)
    return str(value)


def short_summary(value: Any, max_len: int = 85) -> str:
    text = format_value(value, compact=True)
    if len(text) <= max_len:
        return text
    return text[: max_len - 3].rstrip() + "..."


def extra_fields(item: dict[str, Any], defined_fields: set[str]) -> dict[str, Any]:
    aliases = category_aliases()
    extras: dict[str, Any] = {}

    def visit(node: Any) -> None:
        if isinstance(node, dict):
            for key, value in node.items():
                if key in INTERNAL_KEYS or key in aliases:
                    if isinstance(value, (dict, list)):
                        visit(value)
                    continue
                if key not in defined_fields:
                    extras.setdefault(key, value)
                if isinstance(value, (dict, list)):
                    visit(value)
        elif isinstance(node, list):
            for entry in node:
                visit(entry)

    visit(item)
    return extras


def wrap_text(text: str, max_chars: int) -> list[str]:
    words = text.split()
    lines: list[str] = []
    current = ""
    for word in words:
        candidate = word if not current else f"{current} {word}"
        if len(candidate) <= max_chars:
            current = candidate
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)
    return lines or [text]


def wrap_text_preserving_breaks(text: str, max_chars: int) -> list[str]:
    paragraphs = text.splitlines() or [text]
    lines: list[str] = []
    for paragraph in paragraphs:
        stripped = paragraph.strip()
        if not stripped:
            continue
        lines.extend(wrap_text(stripped, max_chars))
    return lines or [text.strip()]


def fit_text_to_box(
    text: str,
    width: float,
    height: float,
    *,
    max_size: int,
    min_size: int = 10,
    padding_x: int = 12,
    padding_y: int = 10,
    line_height_ratio: float = 1.3,
) -> tuple[int, list[str], float]:
    usable_w = max(width - 2 * padding_x, 40)
    usable_h = max(height - 2 * padding_y, 20)
    for size in range(max_size, min_size - 1, -1):
        max_chars = max(4, int(usable_w / max(size * 0.58, 1)))
        lines = wrap_text_preserving_breaks(text, max_chars)
        line_height = max(size * line_height_ratio, size + 2)
        total_height = size + max(0, len(lines) - 1) * line_height
        longest_line = max((len(line) for line in lines), default=0)
        estimated_width = longest_line * size * 0.58
        if total_height <= usable_h and estimated_width <= usable_w:
            return size, lines, line_height
    fallback_size = min_size
    fallback_chars = max(4, int(usable_w / max(fallback_size * 0.58, 1)))
    fallback_lines = wrap_text_preserving_breaks(text, fallback_chars)
    fallback_line_height = max(fallback_size * line_height_ratio, fallback_size + 2)
    return fallback_size, fallback_lines, fallback_line_height


def svg_text(
    x: float,
    y: float,
    text: str,
    *,
    size: int = 16,
    weight: str = "400",
    fill: str = "#0f172a",
    anchor: str = "start",
    max_chars: int | None = None,
    line_height: float | None = None,
) -> str:
    lines = wrap_text(text, max_chars) if max_chars else [text]
    line_height = line_height or size * 1.35
    tspans = []
    for index, line in enumerate(lines):
        dy = "0" if index == 0 else f"{line_height}"
        tspans.append(
            f'<tspan x="{x}" dy="{dy}">{escape(line)}</tspan>'
        )
    return (
        f'<text x="{x}" y="{y}" font-family="Segoe UI, Arial, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" fill="{fill}" text-anchor="{anchor}">'
        + "".join(tspans)
        + "</text>"
    )


def svg_text_box(
    x: float,
    y: float,
    width: float,
    height: float,
    text: str,
    *,
    max_size: int = 16,
    min_size: int = 10,
    weight: str = "400",
    fill: str = "#0f172a",
    align: str = "center",
    valign: str = "middle",
    padding_x: int = 12,
    padding_y: int = 10,
    line_height_ratio: float = 1.3,
) -> str:
    size, lines, line_height = fit_text_to_box(
        text,
        width,
        height,
        max_size=max_size,
        min_size=min_size,
        padding_x=padding_x,
        padding_y=padding_y,
        line_height_ratio=line_height_ratio,
    )
    if align == "center":
        text_x = x + width / 2
        anchor = "middle"
    else:
        text_x = x + padding_x
        anchor = "start"
    total_height = size + max(0, len(lines) - 1) * line_height
    if valign == "middle":
        text_y = y + (height - total_height) / 2 + size * 0.9
    else:
        text_y = y + padding_y + size * 0.9
    tspans = []
    for index, line in enumerate(lines):
        dy = "0" if index == 0 else f"{line_height}"
        tspans.append(f'<tspan x="{text_x}" dy="{dy}">{escape(line)}</tspan>')
    return (
        f'<text x="{text_x}" y="{text_y}" font-family="Segoe UI, Arial, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" fill="{fill}" text-anchor="{anchor}">'
        + "".join(tspans)
        + "</text>"
    )


def svg_rect(
    x: float,
    y: float,
    width: float,
    height: float,
    *,
    fill: str,
    stroke: str = "none",
    stroke_width: int = 1,
    rx: int = 18,
    opacity: float | None = None,
    dash: str | None = None,
) -> str:
    extra = []
    if opacity is not None:
        extra.append(f'opacity="{opacity}"')
    if dash:
        extra.append(f'stroke-dasharray="{dash}"')
    extra_str = " ".join(extra)
    return (
        f'<rect x="{x}" y="{y}" width="{width}" height="{height}" rx="{rx}" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="{stroke_width}" {extra_str}/>'
    )


def svg_line(
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    *,
    stroke: str = "#64748b",
    stroke_width: int = 3,
    dash: str | None = None,
    marker_end: str | None = "url(#arrow)",
) -> str:
    attrs = [
        f'x1="{x1}"',
        f'y1="{y1}"',
        f'x2="{x2}"',
        f'y2="{y2}"',
        f'stroke="{stroke}"',
        f'stroke-width="{stroke_width}"',
        'fill="none"',
    ]
    if dash:
        attrs.append(f'stroke-dasharray="{dash}"')
    if marker_end:
        attrs.append(f'marker-end="{marker_end}"')
    return "<line " + " ".join(attrs) + "/>"


def svg_card(
    x: float,
    y: float,
    width: float,
    height: float,
    title: str,
    body: str,
    *,
    fill: str,
    accent: str,
    title_size: int = 14,
    body_size: int = 22,
) -> str:
    return "\n".join(
        [
            svg_rect(x, y, width, height, fill=fill, stroke=accent, stroke_width=2, rx=22),
            svg_rect(x, y, 8, height, fill=accent, rx=8),
            svg_text_box(x + 18, y + 14, width - 28, 30, title, max_size=title_size, min_size=11, weight="600", fill="#334155"),
            svg_text_box(x + 18, y + 46, width - 28, height - 54, body, max_size=body_size, min_size=14, weight="700", fill="#0f172a"),
        ]
    )


def write_svg(path: Path, width: int, height: int, body: str) -> None:
    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" role="img">'
        "<defs>"
        '<marker id="arrow" markerWidth="10" markerHeight="10" refX="8" refY="5" orient="auto">'
        '<path d="M 0 0 L 10 5 L 0 10 z" fill="#64748b"/>'
        "</marker>"
        "</defs>"
        f"{body}</svg>"
    )
    path.write_text(svg, encoding="utf-8")


def create_assay_heatmap() -> Path:
    path = VISUALS_DIR / "assay_fit_heatmap.svg"
    VISUALS_DIR.mkdir(parents=True, exist_ok=True)
    width, height = 980, 360
    rows = [
        ("cfWGS", [2, 3, 1]),
        ("VirCapSeq on-target", [3, 1, 1]),
        ("VirCapSeq off-target", [1, 1, 1]),
        ("cfRNA / RNA-seq", [1, 1, 3]),
    ]
    cols = [
        "Exogenous virus detection",
        "Endogenous retroelement enumeration",
        "Direct activity readout",
    ]
    colors = {
        1: ("#fee2e2", "#be123c", "Low"),
        2: ("#fef3c7", "#b45309", "Moderate"),
        3: ("#dcfce7", "#166534", "High"),
    }

    body_parts = [
        svg_rect(0, 0, width, height, fill="#f8fafc", rx=0),
        svg_text(40, 42, "Assay Fit by Question", size=28, weight="700"),
        svg_text(
            40,
            72,
            "High-level summary of which assay answers which part of the retrovirus / HERV problem best.",
            size=14,
            fill="#475569",
            max_chars=90,
        ),
    ]

    x0, y0 = 220, 120
    cell_w, cell_h = 210, 52
    for index, col in enumerate(cols):
        x = x0 + index * cell_w
        body_parts.append(svg_text(x + cell_w / 2, y0 - 24, col, size=13, weight="600", fill="#334155", anchor="middle", max_chars=18))

    for row_index, (name, scores) in enumerate(rows):
        y = y0 + row_index * cell_h
        body_parts.append(svg_text(40, y + 34, name, size=15, weight="600", fill="#0f172a", max_chars=20))
        for col_index, score in enumerate(scores):
            x = x0 + col_index * cell_w
            fill, stroke, label = colors[score]
            body_parts.append(svg_rect(x, y, cell_w - 14, cell_h - 10, fill=fill, stroke=stroke, stroke_width=2, rx=16))
            body_parts.append(svg_text(x + (cell_w - 14) / 2, y + 23, label, size=14, weight="700", fill=stroke, anchor="middle"))

    legend_y = y0 + len(rows) * cell_h + 32
    for idx, score in enumerate([3, 2, 1]):
        fill, stroke, label = colors[score]
        lx = 40 + idx * 170
        body_parts.append(svg_rect(lx, legend_y, 28, 28, fill=fill, stroke=stroke, stroke_width=2, rx=8))
        body_parts.append(svg_text(lx + 40, legend_y + 20, label, size=14, fill="#334155"))

    write_svg(path, width, height, "".join(body_parts))
    return path


def create_workflow_diagram() -> Path:
    path = VISUALS_DIR / "analysis_workflow.svg"
    VISUALS_DIR.mkdir(parents=True, exist_ok=True)
    width, height = 1120, 620
    body_parts = [
        svg_rect(0, 0, width, height, fill="#f8fafc", rx=0),
        svg_text(40, 42, "Recommended Analysis Workflow", size=28, weight="700"),
        svg_text(
            40,
            72,
            "The report's practical design: separate exogenous detection, endogenous enumeration, and activity follow-up.",
            size=14,
            fill="#475569",
            max_chars=100,
        ),
    ]

    boxes = {
        "cfWGS": (60, 120, 220, 88, "#dbeafe", "#2563eb", "cfWGS", "Primary substrate for endogenous discovery"),
        "VirCapSeq": (60, 260, 220, 88, "#ede9fe", "#7c3aed", "VirCapSeq", "Best for on-target exogenous detection"),
        "cfRNA": (60, 400, 220, 88, "#dcfce7", "#16a34a", "cfRNA / RNA-seq", "Optional layer for direct activity readout"),
        "Exogenous": (410, 120, 290, 108, "#eff6ff", "#2563eb", "Exogenous virus track", "Competitive alignment plus unique viral or host-virus junction evidence"),
        "Endogenous": (410, 280, 290, 108, "#fefce8", "#ca8a04", "Endogenous retroelement track", "Family-level counts plus xTea and ERVcaller or RetroSnake on deep cfWGS"),
        "Activity": (410, 440, 290, 108, "#f0fdf4", "#16a34a", "Activity follow-up", "Fragmentomics or methylation on cfDNA, or Telescope / ERVmap on RNA"),
        "Outputs": (820, 210, 240, 220, "#fff7ed", "#ea580c", "Decision outputs", "Confident exogenous calls\nCandidate HERV / LINE1 burden\nProxy or direct activity readout"),
    }

    for _, (x, y, w, h, fill, stroke, title, subtitle) in boxes.items():
        body_parts.append(svg_rect(x, y, w, h, fill=fill, stroke=stroke, stroke_width=2, rx=26))
        body_parts.append(svg_text_box(x + 18, y + 14, w - 36, 28, title, max_size=18, min_size=13, weight="700", fill="#0f172a"))
        body_parts.append(svg_text_box(x + 18, y + 42, w - 36, h - 52, subtitle, max_size=13, min_size=10, fill="#334155"))

    body_parts.extend(
        [
            svg_line(280, 164, 410, 164),
            svg_line(280, 304, 410, 174),
            svg_line(280, 304, 410, 334),
            svg_line(280, 444, 410, 494),
            svg_line(700, 174, 820, 250),
            svg_line(700, 334, 820, 320),
            svg_line(700, 494, 820, 390),
            svg_text(312, 272, "off-target host reads only", size=12, fill="#7c3aed", max_chars=20),
            svg_text(312, 470, "optional orthogonal layer", size=12, fill="#16a34a", max_chars=20),
        ]
    )

    write_svg(path, width, height, "".join(body_parts))
    return path


def create_email_answer_grid() -> Path:
    path = VISUALS_DIR / "narrative_email_answer_grid.svg"
    VISUALS_DIR.mkdir(parents=True, exist_ok=True)
    width, height = 1120, 620
    body_parts = [
        svg_rect(0, 0, width, height, fill="#f8fafc", rx=0),
        svg_text(40, 44, "Direct Answer to the Email", size=30, weight="700"),
        svg_text(
            40,
            76,
            "Three decision-focused answers: discrimination, cfDNA adaptation, and assay split.",
            size=14,
            fill="#475569",
            max_chars=98,
        ),
    ]

    cards = [
        (
            40,
            120,
            1040,
            120,
            "#eff6ff",
            "#2563eb",
            "1. Distinguish HIV-1 / HTLV-1 from HERV / LINE1",
            "Yes, but only with a combined host-plus-virus reference, unique viral or host-virus junction evidence, and an explicit ambiguous retroviral bin.",
        ),
        (
            40,
            265,
            1040,
            120,
            "#fefce8",
            "#ca8a04",
            "2. Adapt existing HERV pipelines to cfDNA",
            "Partial transfer is realistic: family-level counting and whitelist-restricted loci are defensible; de novo insertion calling belongs only in deeper cfWGS.",
        ),
        (
            40,
            410,
            1040,
            120,
            "#ecfdf5",
            "#059669",
            "3. Separate WGS from VirCapSeq off-target analysis",
            "Use cfWGS for endogenous HERV and LINE1 work, use VirCapSeq primarily for exogenous HIV-1 and HTLV-1, and treat off-target host reads as secondary exploratory signal only.",
        ),
    ]

    for x, y, w, h, fill, stroke, title, body in cards:
        body_parts.append(svg_rect(x, y, w, h, fill=fill, stroke=stroke, stroke_width=2, rx=28))
        body_parts.append(svg_rect(x, y, 10, h, fill=stroke, rx=10))
        body_parts.append(svg_text_box(x + 24, y + 16, w - 48, 34, title, max_size=20, min_size=13, weight="700", fill="#0f172a"))
        body_parts.append(svg_text_box(x + 24, y + 52, w - 48, h - 64, body, max_size=14, min_size=11, fill="#334155"))

    body_parts.append(
        svg_rect(40, 555, 1040, 40, fill="#fff7ed", stroke="#ea580c", stroke_width=2, rx=16)
    )
    body_parts.append(
        svg_text_box(
            60,
            561,
            1000,
            28,
            "Known HIV-1+ and HTLV-1+ samples are not just validation material; they should serve as calibration controls for the classifier and threshold design.",
            max_size=14,
            min_size=10,
            weight="600",
            fill="#9a3412",
        )
    )

    write_svg(path, width, height, "".join(body_parts))
    return path


def create_discrimination_logic_diagram() -> Path:
    path = VISUALS_DIR / "narrative_discrimination_logic.svg"
    VISUALS_DIR.mkdir(parents=True, exist_ok=True)
    width, height = 1180, 620
    body_parts = [
        svg_rect(0, 0, width, height, fill="#f8fafc", rx=0),
        svg_text(40, 44, "Discrimination Logic for Retrovirus-like Reads", size=30, weight="700"),
        svg_text(
            40,
            76,
            "The key design choice is to preserve ambiguity rather than forcing every read into a viral or endogenous bin.",
            size=14,
            fill="#475569",
            max_chars=104,
        ),
    ]

    body_parts.extend(
        [
            svg_rect(50, 170, 260, 100, fill="#e0f2fe", stroke="#0284c7", stroke_width=2, rx=24),
            svg_text_box(70, 182, 220, 26, "Input reads", max_size=20, min_size=14, weight="700", fill="#0f172a"),
            svg_text_box(70, 214, 220, 42, "cfWGS or VirCapSeq retrovirus-like fragments", max_size=14, min_size=11, fill="#334155"),
            svg_rect(360, 130, 330, 180, fill="#eef2ff", stroke="#4f46e5", stroke_width=2, rx=24),
            svg_text_box(382, 144, 286, 28, "Combined reference", max_size=20, min_size=14, weight="700", fill="#0f172a"),
            svg_text_box(382, 178, 286, 42, "hg38 + HIV-1 + HTLV-1 + HERV / LINE1 decoys", max_size=14, min_size=11, fill="#334155"),
            svg_text_box(382, 228, 286, 56, "Goal: prevent ambiguous reads from being forced into an exogenous call.", max_size=13, min_size=10, fill="#4338ca"),
            svg_rect(740, 130, 390, 180, fill="#f8fafc", stroke="#64748b", stroke_width=2, rx=24),
            svg_text_box(762, 144, 346, 26, "Decision rules", max_size=20, min_size=14, weight="700", fill="#0f172a"),
            svg_text_box(762, 180, 346, 34, "1. Unique viral region or host-virus junction", max_size=14, min_size=11, fill="#334155"),
            svg_text_box(762, 220, 346, 34, "2. Human-flanked endogenous locus support", max_size=14, min_size=11, fill="#334155"),
            svg_text_box(762, 260, 346, 34, "3. Conserved LTR / gag / pol only -> ambiguous", max_size=14, min_size=11, fill="#334155"),
            svg_rect(80, 400, 280, 120, fill="#ecfdf5", stroke="#059669", stroke_width=2, rx=26),
            svg_text_box(102, 414, 236, 28, "Unique exogenous", max_size=22, min_size=14, weight="700", fill="#065f46"),
            svg_text_box(102, 450, 236, 48, "Use for HIV-1 or HTLV-1 calls and proviral follow-up.", max_size=14, min_size=11, fill="#334155"),
            svg_rect(450, 400, 280, 120, fill="#fff7ed", stroke="#ea580c", stroke_width=2, rx=26),
            svg_text_box(472, 414, 236, 28, "Ambiguous retroviral", max_size=22, min_size=14, weight="700", fill="#9a3412"),
            svg_text_box(472, 448, 236, 52, "Keep separate; do not overinterpret as either viral infection or HERV activation.", max_size=14, min_size=11, fill="#334155"),
            svg_rect(820, 400, 280, 120, fill="#eff6ff", stroke="#2563eb", stroke_width=2, rx=26),
            svg_text_box(842, 414, 236, 28, "Unique endogenous", max_size=22, min_size=14, weight="700", fill="#1d4ed8"),
            svg_text_box(842, 450, 236, 48, "Use for HERV or LINE1 burden and selective insertion analyses.", max_size=14, min_size=11, fill="#334155"),
            svg_line(310, 220, 360, 220),
            svg_line(690, 220, 740, 220),
            svg_line(860, 310, 860, 360, stroke="#64748b"),
            svg_line(860, 360, 220, 360, stroke="#64748b", marker_end=None),
            svg_line(220, 360, 220, 400, stroke="#64748b"),
            svg_line(860, 360, 590, 360, stroke="#64748b", marker_end=None),
            svg_line(590, 360, 590, 400, stroke="#64748b"),
            svg_line(860, 360, 960, 360, stroke="#64748b", marker_end=None),
            svg_line(960, 360, 960, 400, stroke="#64748b"),
        ]
    )

    write_svg(path, width, height, "".join(body_parts))
    return path


def create_control_calibration_diagram() -> Path:
    path = VISUALS_DIR / "narrative_control_calibration.svg"
    VISUALS_DIR.mkdir(parents=True, exist_ok=True)
    width, height = 1180, 640
    body_parts = [
        svg_rect(0, 0, width, height, fill="#f8fafc", rx=0),
        svg_text(40, 44, "How the Known Infected Samples Should Be Used", size=30, weight="700"),
        svg_text(
            40,
            76,
            "The positive controls are not peripheral; they define the calibration strategy for the entire pipeline.",
            size=14,
            fill="#475569",
            max_chars=106,
        ),
    ]

    sample_cards = [
        (60, 150, "#fee2e2", "#dc2626", "HIV-1+ plasma", "Defines unique HIV-supporting regions and expected exogenous alignment patterns."),
        (60, 290, "#ffedd5", "#ea580c", "HTLV-1+ plasma", "Defines unique HTLV-supporting regions and integration-oriented evidence when present."),
        (60, 430, "#e0f2fe", "#0284c7", "Negative plasma", "Quantifies background assignment into viral and endogenous bins under the same preprocessing."),
    ]

    for x, y, fill, stroke, title, body in sample_cards:
        body_parts.append(svg_rect(x, y, 270, 120, fill=fill, stroke=stroke, stroke_width=2, rx=24))
        body_parts.append(svg_text_box(x + 16, y + 14, 238, 30, title, max_size=20, min_size=13, weight="700", fill="#0f172a"))
        body_parts.append(svg_text_box(x + 16, y + 46, 238, 58, body, max_size=13, min_size=10, fill="#334155"))

    body_parts.append(svg_rect(420, 230, 320, 230, fill="#eef2ff", stroke="#4f46e5", stroke_width=2, rx=28))
    body_parts.append(svg_text_box(440, 246, 280, 34, "Calibration cohort", max_size=24, min_size=14, weight="700", fill="#312e81"))
    body_parts.append(svg_text_box(440, 290, 280, 76, "Run all controls through the same trimming, alignment, and deduplication settings planned for the pilot.", max_size=14, min_size=10, fill="#334155"))
    body_parts.append(svg_text_box(440, 372, 280, 76, "Measure leakage into HERV / LINE1 decoys and quantify the ambiguous retroviral bin before interpreting real samples.", max_size=14, min_size=10, fill="#334155"))

    output_cards = [
        (850, 140, "#ecfdf5", "#059669", "Discriminative regions", "Which HIV-1 or HTLV-1 segments remain uniquely assignable after short-fragment filtering."),
        (850, 290, "#fefce8", "#ca8a04", "Misassignment estimates", "How often viral reads bleed into HERV / LINE1 decoys and how much background persists in negatives."),
        (850, 440, "#eff6ff", "#2563eb", "Locked decision thresholds", "Rules for unique exogenous, unique endogenous, and ambiguous bins before pilot interpretation."),
    ]
    for x, y, fill, stroke, title, body in output_cards:
        body_parts.append(svg_rect(x, y, 270, 140, fill=fill, stroke=stroke, stroke_width=2, rx=24))
        body_parts.append(svg_text_box(x + 16, y + 16, 238, 42, title, max_size=20, min_size=13, weight="700", fill="#0f172a"))
        body_parts.append(svg_text_box(x + 16, y + 62, 238, 62, body, max_size=13, min_size=10, fill="#334155"))

    body_parts.extend(
        [
            svg_line(330, 210, 420, 305),
            svg_line(330, 350, 420, 345),
            svg_line(330, 490, 420, 385),
            svg_line(740, 305, 850, 210),
            svg_line(740, 345, 850, 360),
            svg_line(740, 385, 850, 510),
        ]
    )

    write_svg(path, width, height, "".join(body_parts))
    return path


def create_narrative_visual_assets() -> None:
    create_discrimination_logic_diagram()
    create_control_calibration_diagram()


def parse_claim_counts() -> Counter[str]:
    counts: Counter[str] = Counter()
    if not CLAIM_VALIDATION_PATH.exists():
        return counts
    for line in CLAIM_VALIDATION_PATH.read_text(encoding="utf-8").splitlines():
        if not line.startswith("| C"):
            continue
        parts = [part.strip() for part in line.split("|")[1:-1]]
        if len(parts) < 3:
            continue
        verdict = parts[2]
        verdict_key = "Direct" if verdict.startswith("Direct") else verdict
        counts[verdict_key] += 1
    return counts


def parse_warning_counts() -> Counter[str]:
    counts: Counter[str] = Counter()
    if not SOURCE_VALIDATION_JSON_PATH.exists():
        return counts
    data = json.loads(SOURCE_VALIDATION_JSON_PATH.read_text(encoding="utf-8"))
    for record in data:
        warning = record.get("warning")
        if warning:
            counts[warning] += 1
    return counts


def create_validation_dashboard() -> Path:
    path = VISUALS_DIR / "validation_dashboard.svg"
    VISUALS_DIR.mkdir(parents=True, exist_ok=True)
    width, height = 1120, 560
    summary = load_yaml(SOURCE_SUMMARY_PATH) if SOURCE_SUMMARY_PATH.exists() else {}
    total_sources = int(summary.get("total_sources", 0))
    warnings = int(summary.get("sources_with_warnings", 0))
    clean = max(total_sources - warnings, 0)
    claim_counts = parse_claim_counts()
    warning_counts = parse_warning_counts()

    body_parts = [
        svg_rect(0, 0, width, height, fill="#f8fafc", rx=0),
        svg_text(40, 42, "Validation Dashboard", size=28, weight="700"),
        svg_text(
            40,
            72,
            "Quick visual readout of the source audit and the strength of the report's main claims.",
            size=14,
            fill="#475569",
            max_chars=96,
        ),
    ]

    cards = [
        ("Sources audited", str(total_sources), "#eff6ff", "#2563eb"),
        ("Clean sources", str(clean), "#f0fdf4", "#16a34a"),
        ("Sources with warnings", str(warnings), "#fff7ed", "#ea580c"),
        ("Claims reviewed", str(sum(claim_counts.values())), "#f5f3ff", "#7c3aed"),
    ]
    card_x = 40
    for title, value, fill, accent in cards:
        body_parts.append(svg_card(card_x, 110, 245, 92, title, value, fill=fill, accent=accent))
        card_x += 265

    body_parts.append(svg_text(40, 250, "Claim verdicts", size=20, weight="700"))
    verdict_colors = {
        "Direct": "#2563eb",
        "Supported + inference": "#7c3aed",
        "Search-based absence claim": "#ea580c",
    }
    max_claim = max(claim_counts.values()) if claim_counts else 1
    start_y = 290
    for idx, verdict in enumerate(["Direct", "Supported + inference", "Search-based absence claim"]):
        count = claim_counts.get(verdict, 0)
        y = start_y + idx * 62
        bar_w = 320 * count / max_claim if max_claim else 0
        body_parts.append(svg_text(40, y + 18, verdict, size=14, weight="600", fill="#334155", max_chars=24))
        body_parts.append(svg_rect(240, y, 340, 28, fill="#e2e8f0", rx=14))
        body_parts.append(svg_rect(240, y, bar_w, 28, fill=verdict_colors[verdict], rx=14))
        body_parts.append(svg_text(595, y + 20, str(count), size=14, weight="700", fill="#0f172a"))

    body_parts.append(svg_text(640, 250, "Warning types", size=20, weight="700"))
    max_warning = max(warning_counts.values()) if warning_counts else 1
    for idx, warning in enumerate(sorted(warning_counts)):
        count = warning_counts[warning]
        y = 290 + idx * 62
        label = warning.replace("_", " ")
        bar_w = 320 * count / max_warning if max_warning else 0
        body_parts.append(svg_text(640, y + 18, label, size=14, weight="600", fill="#334155", max_chars=24))
        body_parts.append(svg_rect(860, y, 180, 28, fill="#e2e8f0", rx=14))
        body_parts.append(svg_rect(860, y, bar_w, 28, fill="#ef4444" if "http" in warning else "#f59e0b", rx=14))
        body_parts.append(svg_text(1050, y + 20, str(count), size=14, weight="700", fill="#0f172a"))

    body_parts.append(
        svg_text(
            40,
            515,
            "Interpretation: most remaining source warnings are access-control artifacts rather than contradictions in the underlying tool or paper metadata.",
            size=13,
            fill="#475569",
            max_chars=120,
        )
    )

    write_svg(path, width, height, "".join(body_parts))
    return path


def render_visual_overview() -> list[str]:
    heatmap = create_assay_heatmap()
    workflow = create_workflow_diagram()
    dashboard = create_validation_dashboard()
    return [
        "## Visual Overview",
        "",
        "### Assay Fit by Question",
        "",
        local_image(heatmap, "Assay fit heatmap"),
        "",
        "### Recommended Workflow",
        "",
        local_image(workflow, "Analysis workflow"),
        "",
        "### Validation Dashboard",
        "",
        local_image(dashboard, "Validation dashboard"),
        "",
    ]


def extract_docx_paragraphs(path: Path) -> list[str]:
    ns = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
    with zipfile.ZipFile(path) as zf:
        xml = zf.read("word/document.xml")
    root = ET.fromstring(xml)
    paragraphs = []
    for para in root.findall(".//w:p", ns):
        runs = [t.text for t in para.findall(".//w:t", ns) if t.text]
        if runs:
            paragraphs.append("".join(runs).strip())
    return [p for p in paragraphs if p]


def parse_percent(value: str) -> float | None:
    value = value.strip()
    if not value or value.upper() == "NA":
        return None
    if ":" in value:
        value = value.split(":", 1)[1].strip()
    value = value.rstrip("%")
    try:
        return float(value)
    except ValueError:
        return None


def parse_cpu_hours(value: str) -> float | None:
    value = value.strip()
    match = re.fullmatch(r"(\d+):(\d+):(\d+)", value)
    if not match:
        return None
    hours, minutes, seconds = (int(part) for part in match.groups())
    return hours + minutes / 60 + seconds / 3600


def parse_int_strict(value: str) -> int | None:
    value = value.strip()
    if not value or value.upper() == "NA":
        return None
    return int(value) if re.fullmatch(r"\d+", value) else None


def load_benchmark_table() -> list[dict[str, Any]]:
    if not BENCHMARK_DOCX_PATH.exists():
        return []
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    paragraphs = extract_docx_paragraphs(BENCHMARK_DOCX_PATH)
    rows = []
    tools = ["Retroseq", "Melt", "Retroseq+", "Steak", "ERVcaller"]
    tool_classes = {
        "Retroseq": "Generalist",
        "Melt": "Generalist",
        "Retroseq+": "HERV-specific",
        "Steak": "HERV-specific",
        "ERVcaller": "HERV-specific",
    }
    try:
        start = paragraphs.index("Retroseq")
    except ValueError:
        return []

    cursor = start
    for tool in tools:
        if cursor >= len(paragraphs) or paragraphs[cursor] != tool:
            break
        long_read_text = paragraphs[cursor + 4]
        offset = 0
        if cursor + 5 < len(paragraphs):
            next_value = paragraphs[cursor + 5]
            if parse_cpu_hours(next_value) is None and "%" in next_value:
                long_read_text = f"{long_read_text}; {next_value}"
                offset = 1
        rows.append(
            {
                "tool": tool,
                "tool_class": tool_classes[tool],
                "total_predictions": int(paragraphs[cursor + 1]),
                "avg_sensitivity": float(paragraphs[cursor + 2]),
                "avg_precision": float(paragraphs[cursor + 3]),
                "long_read_confirmed_text": long_read_text,
                "long_read_confirmed_percent": parse_percent(long_read_text),
                "cpu_time_text": paragraphs[cursor + 5 + offset],
                "cpu_time_hours": parse_cpu_hours(paragraphs[cursor + 5 + offset]),
                "citations": parse_int_strict(paragraphs[cursor + 6 + offset]),
                "year_published": parse_int_strict(paragraphs[cursor + 7 + offset]),
            }
        )
        cursor += 8 + offset

    csv_lines = [
        "tool,tool_class,total_predictions,avg_sensitivity,avg_precision,long_read_confirmed_text,long_read_confirmed_percent,cpu_time_text,cpu_time_hours,citations,year_published"
    ]
    for row in rows:
        csv_lines.append(
            ",".join(
                [
                    row["tool"],
                    row["tool_class"],
                    str(row["total_predictions"]),
                    str(row["avg_sensitivity"]),
                    str(row["avg_precision"]),
                    '"' + row["long_read_confirmed_text"].replace('"', '""') + '"',
                    "" if row["long_read_confirmed_percent"] is None else str(row["long_read_confirmed_percent"]),
                    row["cpu_time_text"],
                    "" if row["cpu_time_hours"] is None else f"{row['cpu_time_hours']:.3f}",
                    "" if row["citations"] is None else str(row["citations"]),
                    "" if row["year_published"] is None else str(row["year_published"]),
                ]
            )
        )
    BENCHMARK_CSV_PATH.write_text("\n".join(csv_lines) + "\n", encoding="utf-8")
    return rows


def create_benchmark_plot(rows: list[dict[str, Any]]) -> Path | None:
    if not rows:
        return None
    path = VISUALS_DIR / "herv_tool_benchmark.svg"
    VISUALS_DIR.mkdir(parents=True, exist_ok=True)
    width, height = 1280, 700
    colors = {"HERV-specific": "#2563eb", "Generalist": "#f97316"}
    metric_bg = "#e2e8f0"
    sens_x = 300
    prec_x = 630
    bar_w = 220

    body_parts = [
        svg_rect(0, 0, width, height, fill="#f8fafc", rx=0),
        svg_text(40, 42, "Benchmark Data from Supplementary Table 3", size=28, weight="700"),
        svg_text(
            40,
            72,
            "Downloaded from the benchmark repository associated with the 2022 tool-assessment paper.",
            size=14,
            fill="#475569",
            max_chars=96,
        ),
        svg_text(40, 652, f"Source data: {BENCHMARK_CSV_PATH.name}", size=13, fill="#64748b"),
    ]

    body_parts.extend(
        [
            svg_text(300, 118, "Sensitivity", size=18, weight="700", fill="#0f172a"),
            svg_text(630, 118, "Precision", size=18, weight="700", fill="#0f172a"),
            svg_text(980, 118, "Interpretation", size=22, weight="700", fill="#0f172a"),
        ]
    )

    for tick in [0.0, 0.25, 0.5, 0.75, 1.0]:
        sx = sens_x + tick * bar_w
        px = prec_x + tick * bar_w
        body_parts.append(svg_line(sx, 140, sx, 540, marker_end=None, stroke="#e2e8f0", stroke_width=1))
        body_parts.append(svg_line(px, 140, px, 540, marker_end=None, stroke="#e2e8f0", stroke_width=1))
        body_parts.append(svg_text(sx, 560, f"{tick:.2f}".rstrip("0").rstrip("."), size=12, fill="#64748b", anchor="middle"))
        body_parts.append(svg_text(px, 560, f"{tick:.2f}".rstrip("0").rstrip("."), size=12, fill="#64748b", anchor="middle"))

    row_y = 170
    row_gap = 72
    sorted_rows = sorted(rows, key=lambda row: (-row["avg_sensitivity"], -row["avg_precision"], row["tool"]))
    for idx, row in enumerate(sorted_rows):
        y = row_y + idx * row_gap
        color = colors[row["tool_class"]]
        body_parts.append(svg_text(40, y + 20, row["tool"], size=18, weight="700", fill="#0f172a", max_chars=18))
        body_parts.append(svg_rect(40, y + 30, 124, 26, fill="#ffffff", stroke=color, stroke_width=2, rx=10))
        body_parts.append(svg_text(52, y + 48, row["tool_class"], size=12, weight="600", fill=color, max_chars=18))
        body_parts.append(svg_text(180, y + 48, f"{row['total_predictions']} predictions", size=12, fill="#64748b"))

        body_parts.append(svg_rect(sens_x, y, bar_w, 28, fill=metric_bg, rx=14))
        body_parts.append(svg_rect(sens_x, y, bar_w * row["avg_sensitivity"], 28, fill="#2563eb", rx=14, opacity=0.85))
        body_parts.append(svg_text(sens_x + bar_w + 16, y + 20, f"{row['avg_sensitivity']:.2f}", size=13, weight="700", fill="#1e3a8a"))

        body_parts.append(svg_rect(prec_x, y, bar_w, 28, fill=metric_bg, rx=14))
        body_parts.append(svg_rect(prec_x, y, bar_w * row["avg_precision"], 28, fill="#0f766e" if row["tool"] == "ERVcaller" else "#334155", rx=14, opacity=0.85))
        body_parts.append(svg_text(prec_x + bar_w + 16, y + 20, f"{row['avg_precision']:.2f}", size=13, weight="700", fill="#0f172a"))

    best_sens = max(rows, key=lambda row: row["avg_sensitivity"])
    best_prec = max(rows, key=lambda row: row["avg_precision"])
    best_confirmed = max(
        [row for row in rows if row["long_read_confirmed_percent"] is not None],
        key=lambda row: row["long_read_confirmed_percent"],
    )

    body_parts.append(svg_rect(960, 150, 260, 132, fill="#ffffff", stroke="#2563eb", stroke_width=2, rx=24))
    body_parts.append(svg_text_box(978, 162, 224, 28, "Tool-specific takeaways", max_size=18, min_size=12, weight="700", fill="#0f172a"))
    body_parts.append(svg_text_box(978, 194, 224, 72, f"Highest sensitivity: {best_sens['tool']} ({best_sens['avg_sensitivity']:.2f})\nHighest precision: {best_prec['tool']} ({best_prec['avg_precision']:.2f})", max_size=14, min_size=10, fill="#334155"))

    body_parts.append(svg_rect(960, 308, 260, 120, fill="#ffffff", stroke="#ea580c", stroke_width=2, rx=24))
    body_parts.append(svg_text_box(978, 320, 224, 28, "Confirmation note", max_size=18, min_size=12, weight="700", fill="#0f172a"))
    body_parts.append(svg_text_box(978, 352, 224, 56, f"Highest single-value long-read confirmation: {best_confirmed['tool']} ({best_confirmed['long_read_confirmed_percent']:.0f}%)", max_size=14, min_size=10, fill="#334155"))

    body_parts.append(svg_rect(960, 452, 260, 118, fill="#ffffff", stroke="#64748b", stroke_width=2, rx=24))
    body_parts.append(svg_text_box(978, 462, 224, 28, "What this supports", max_size=18, min_size=12, weight="700", fill="#0f172a"))
    body_parts.append(svg_text_box(978, 494, 224, 60, "Supports ERVcaller plus companion-caller confirmation, not a blanket HERV-specific-over-generalist rule.", max_size=12, min_size=9, fill="#334155"))

    body_parts.append(svg_text(300, 592, "Average sensitivity on simulated HML-6 insertions", size=13, weight="600", fill="#475569", max_chars=34))
    body_parts.append(svg_text(630, 592, "Average precision on simulated HML-6 insertions", size=13, weight="600", fill="#475569", max_chars=34))

    write_svg(path, width, height, "".join(body_parts))
    return path


def render_empirical_section(rows: list[dict[str, Any]]) -> list[str]:
    lines = ["## Downloaded Benchmark Support", ""]
    if not rows:
        lines.append("No benchmark table could be extracted from the downloaded supplementary repository.")
        lines.append("")
        return lines

    plot_path = create_benchmark_plot(rows)
    lines.append(
        "I downloaded the benchmark repository linked to the HERV tool-assessment paper and extracted Supplementary Table 3 into "
        + local_md_link(BENCHMARK_CSV_PATH)
        + "."
    )
    lines.append("")
    lines.append("This benchmark supports the report's tool-specific ranking, especially the decision to emphasize `ERVcaller` and to treat `STEAK` as a secondary comparator rather than a sole discovery tool.")
    lines.append("")
    if plot_path:
        lines.append(local_image(plot_path, "Downloaded HERV benchmark plot"))
        lines.append("")

    herv_specific = [row for row in rows if row["tool_class"] == "HERV-specific"]
    generalist = [row for row in rows if row["tool_class"] == "Generalist"]
    mean_hs = sum(row["avg_sensitivity"] for row in herv_specific) / len(herv_specific)
    mean_hp = sum(row["avg_precision"] for row in herv_specific) / len(herv_specific)
    mean_gs = sum(row["avg_sensitivity"] for row in generalist) / len(generalist)
    mean_gp = sum(row["avg_precision"] for row in generalist) / len(generalist)
    best_sens = max(rows, key=lambda row: row["avg_sensitivity"])
    best_prec = max(rows, key=lambda row: row["avg_precision"])
    best_confirmed = max(
        [row for row in rows if row["long_read_confirmed_percent"] is not None],
        key=lambda row: row["long_read_confirmed_percent"],
    )

    lines.append(f"- Highest sensitivity in the downloaded table: `{best_sens['tool']}` at `{best_sens['avg_sensitivity']:.2f}`")
    lines.append(f"- Highest precision in the downloaded table: `{best_prec['tool']}` at `{best_prec['avg_precision']:.2f}`")
    lines.append(f"- Highest long-read confirmation among single-value rows: `{best_confirmed['tool']}` at `{best_confirmed['long_read_confirmed_percent']:.0f}%`")
    lines.append(f"- Class averages are mixed: HERV-specific mean sensitivity `{mean_hs:.2f}` vs generalist `{mean_gs:.2f}`, largely because `STEAK` is precision-heavy but sensitivity-poor")
    lines.append(f"- HERV-specific mean precision remains slightly higher: `{mean_hp:.2f}` vs generalist `{mean_gp:.2f}`")
    lines.append("")
    lines.append("The strongest support here is tool-specific rather than class-average: `ERVcaller` has the best balance in this benchmark, `STEAK` fits the report's characterization as a secondary high-precision comparator, and `Retroseq+` contributes confirmation-style value rather than leading on sensitivity or precision.")
    lines.append("")
    lines.append("This downloaded benchmark does not by itself establish cfDNA suitability, which remains an inferential extension from WGS-focused benchmarking plus cfDNA assay constraints.")
    lines.append("")
    return lines


def render_email_aligned_summary(rows: list[dict[str, Any]]) -> list[str]:
    lines = ["## Direct Answer to the Email", ""]
    lines.append("This report is organized around the three concrete questions in the email, with the item-by-item deep-research summaries moved below as supporting material.")
    lines.append("")
    lines.append("### 1. Can HIV-1 and HTLV-1 be distinguished from HERVs and LINE1?")
    lines.append("")
    lines.append("Yes, but only with an explicit ambiguity-aware discrimination framework. The recommended rule is to align against a combined reference containing `hg38`, HIV-1, HTLV-1, and endogenous `HERV/LINE1` decoys, call `exogenous` only from unique viral regions or host-virus junctions, call `endogenous` only from human-flanked loci, and keep a third `ambiguous retroviral` bin for conserved retroviral sequence.")
    lines.append("")
    lines.append("### 2. Can existing HERV enumeration pipelines be adapted to cfDNA?")
    lines.append("")
    lines.append("Yes, but only partially. Existing HERV DNA callers can support selective insertion calling from sufficiently deep `cfWGS`, but the safest cfDNA transfer is family-level counting plus whitelist-restricted locus analysis. They should not be described as validated drop-in cfDNA pipelines.")
    lines.append("")
    lines.append("### 3. How should this be handled for WGS versus VirCapSeq off-target reads?")
    lines.append("")
    lines.append("Treat them as different substrates. `cfWGS` should be the primary dataset for endogenous `HERV/LINE1` work. `VirCapSeq` should be used primarily for exogenous HIV-1 and HTLV-1 detection, while off-target host reads should be limited to secondary exploratory summaries rather than used as the main discovery substrate for new endogenous insertions.")
    lines.append("")
    lines.append("## Recommended Pilot Plan")
    lines.append("")
    lines.append("1. Use the known `HIV-1+`, `HTLV-1+`, and negative samples as a calibration cohort, not just as retrospective validation material.")
    lines.append("2. Build a combined host-plus-virus reference and predefine the three-bin classification rule: `unique exogenous`, `unique endogenous`, and `ambiguous retroviral`.")
    lines.append("3. Run the exogenous discrimination track first on all available `cfWGS` and `VirCapSeq` data to quantify viral-versus-endogenous misassignment under the real assay conditions.")
    lines.append("4. Run endogenous enumeration primarily on `cfWGS`, using family-level repeat counting plus `xTea` and a HERV-focused companion caller such as `ERVcaller` or `RetroSnake` on the deeper samples.")
    lines.append("5. Use `VirCapSeq` off-target host reads only if deduplicated host yield is adequate for coarse endogenous summaries; otherwise stop at the exogenous-retrovirus result.")
    lines.append("")
    if rows:
        best_sens = max(rows, key=lambda row: row["avg_sensitivity"])
        lines.append(
            "The downloaded benchmark supports the tool-specific emphasis in this plan: "
            f"`{best_sens['tool']}` has the highest sensitivity in the extracted table "
            f"at `{best_sens['avg_sensitivity']:.2f}`, which is why the report favors `ERVcaller`-style sensitivity plus companion-caller confirmation over a single rigid caller choice."
        )
        lines.append("")
    lines.append("`cfRNA`, methylation, and fragmentomics remain reasonable follow-on extensions if the question later shifts from enumeration to activity, but they are secondary to the present email's main request.")
    lines.append("")
    return lines


def render_toc(items: list[dict[str, Any]]) -> list[str]:
    lines = ["## Table of Contents", ""]
    for index, item in enumerate(items, start=1):
        anchor = slugify(str(item["name"]))
        summary_parts = []
        for field in SUMMARY_FIELDS:
            value = item["fields"].get(field)
            if field in item["uncertain"] or is_uncertain_value(value):
                continue
            label = field.replace("_", " ").title()
            summary_parts.append(f"{label}: {short_summary(value)}")
        summary = f" - {' | '.join(summary_parts)}" if summary_parts else ""
        lines.append(f"{index}. [{item['name']}](#{anchor}){summary}")
    lines.append("")
    return lines


def render_validation_snapshot() -> list[str]:
    lines = ["## Validation Snapshot", ""]
    summary_path = PROJECT_DIR / "validation" / "source_summary.yaml"
    claim_path = PROJECT_DIR / "validation" / "claim_validation.md"
    source_path = PROJECT_DIR / "validation" / "source_validation.md"
    narrative_path = PROJECT_DIR / "report_narrative.md"

    if summary_path.exists():
        summary = load_yaml(summary_path)
        lines.append(f"- Sources audited: {summary.get('total_sources', 'n/a')}")
        lines.append(f"- HTTP OK: {summary.get('ok_sources', 'n/a')}")
        lines.append(f"- Sources with warnings: {summary.get('sources_with_warnings', 'n/a')}")
    lines.append(f"- Narrative report: {local_md_link(narrative_path)}")
    if claim_path.exists():
        lines.append(f"- Claim validation: {local_md_link(claim_path)}")
    if source_path.exists():
        lines.append(f"- Source audit: {local_md_link(source_path)}")
    lines.append("")
    return lines


def render_item(item: dict[str, Any], field_categories: list[dict[str, Any]]) -> list[str]:
    lines = [f"## {item['name']}", ""]
    for category in field_categories:
        category_name = category["category"]
        category_lines: list[str] = []
        for field in category.get("fields", []):
            name = field["name"]
            if name in item["uncertain"]:
                continue
            value = item["fields"].get(name)
            if is_uncertain_value(value):
                continue
            formatted = format_value(value)
            if not formatted:
                continue
            label = name.replace("_", " ").title()
            if "\n" in formatted:
                category_lines.append(f"### {label}")
                category_lines.append("")
                category_lines.append(formatted)
                category_lines.append("")
            else:
                category_lines.append(f"- **{label}:** {formatted}")
        if category_lines:
            lines.append(f"### {category_name}")
            lines.append("")
            lines.extend(category_lines)
            if lines[-1] != "":
                lines.append("")

    extras = extra_fields(item["raw"], item["defined_fields"])
    filtered_extras = {
        key: value
        for key, value in extras.items()
        if key not in item["uncertain"] and not is_uncertain_value(value)
    }
    if filtered_extras:
        lines.append("### Other Info")
        lines.append("")
        for key, value in filtered_extras.items():
            formatted = format_value(value)
            if not formatted:
                continue
            label = key.replace("_", " ").title()
            if "\n" in formatted:
                lines.append(f"#### {label}")
                lines.append("")
                lines.append(formatted)
                lines.append("")
            else:
                lines.append(f"- **{label}:** {formatted}")
        if lines[-1] != "":
            lines.append("")
    return lines


def load_results(project_dir: Path, results_dir: Path, field_names: set[str]) -> list[dict[str, Any]]:
    items = []
    for json_path in sorted(results_dir.glob("*.json")):
        raw = json.loads(json_path.read_text(encoding="utf-8"))
        if not isinstance(raw, dict):
            continue
        fields = collect_all_fields(raw, field_names)
        uncertain = set(raw.get("uncertain", []))
        name = fields.get("name", json_path.stem)
        items.append(
            {
                "name": name,
                "raw": raw,
                "fields": fields,
                "uncertain": uncertain,
                "defined_fields": field_names,
                "path": json_path,
            }
        )
    return items


def main() -> None:
    outline = load_yaml(OUTLINE_PATH)
    fields_yaml = load_yaml(FIELDS_PATH)
    results_dir = (PROJECT_DIR / outline["execution"]["output_dir"]).resolve()
    field_categories = fields_yaml.get("field_categories", [])
    field_names = {
        field["name"]
        for category in field_categories
        for field in category.get("fields", [])
    }
    items = load_results(PROJECT_DIR, results_dir, field_names)
    benchmark_rows = load_benchmark_table()
    create_narrative_visual_assets()

    lines: list[str] = []
    lines.append(f"# {outline['topic']}")
    lines.append("")
    lines.append("Generated from the deep-research JSON outputs using the `research-report` workflow and revised to match the email's decision-focused framing.")
    lines.append("")
    lines.extend(render_email_aligned_summary(benchmark_rows))
    lines.extend(render_validation_snapshot())
    lines.extend(render_visual_overview())
    lines.extend(render_empirical_section(benchmark_rows))
    lines.extend(render_toc(items))

    for item in items:
        lines.extend(render_item(item, field_categories))

    OUTPUT_PATH.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
