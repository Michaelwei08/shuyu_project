from __future__ import annotations

import argparse
import csv
import re
import subprocess
from pathlib import Path


DEFAULT_BASE_CATEGORIES = ["HUMAN", "HERV", "LINE1"]


def parse_header_id(header: str) -> str:
    return header[1:].strip().split()[0]


def safe_id(text: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", text.strip())
    return cleaned[:80] or "record"


def iter_fasta(path: Path):
    header: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, chunks
                header = line
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, chunks


def load_reference_map(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or "reference_id" not in rows[0] or "category" not in rows[0]:
        raise ValueError("base reference map must contain reference_id and category")
    return {row["reference_id"]: row for row in rows}


def infer_panel_category(header: str) -> str:
    text = re.sub(r"[^A-Z0-9]+", " ", header.upper())
    compact = text.replace(" ", "")
    if "HIV1" in compact or "HIV 1" in text or "IMMUNODEFICIENCY VIRUS 1" in text:
        return "HIV1"
    if "HIV2" in compact or "HIV 2" in text or "IMMUNODEFICIENCY VIRUS 2" in text:
        return "HIV2"
    if "HTLV1" in compact or "HTLV I" in text or "LYMPHOTROPIC VIRUS 1" in text:
        return "HTLV1"
    if "HTLV2" in compact or "HTLV II" in text or "LYMPHOTROPIC VIRUS 2" in text:
        return "HTLV2"
    return "OTHER_VIRAL"


def write_record(handle, header: str, chunks: list[str]) -> None:
    handle.write(header + "\n")
    for chunk in chunks:
        if chunk:
            handle.write(chunk.upper() + "\n")


def append_base_records(
    base_fasta: Path,
    base_map: dict[str, dict[str, str]],
    categories: set[str],
    output_handle,
    map_rows: list[dict[str, str]],
) -> None:
    written = 0
    for header, chunks in iter_fasta(base_fasta):
        ref_id = parse_header_id(header)
        row = base_map.get(ref_id)
        if row is None or row["category"] not in categories:
            continue
        write_record(output_handle, header, chunks)
        map_rows.append(
            {
                "reference_id": ref_id,
                "category": row["category"],
                "source": row.get("source", f"base:{base_fasta}"),
                "description": row.get("description", header[1:].strip()),
                "original_reference_id": ref_id,
                "original_header": header[1:].strip(),
            }
        )
        written += 1
    if written == 0:
        raise ValueError(f"no base records matched categories: {', '.join(sorted(categories))}")


def append_panel_records(
    panel_fasta: Path,
    output_handle,
    map_rows: list[dict[str, str]],
    inventory_rows: list[dict[str, str]],
) -> None:
    for idx, (header, chunks) in enumerate(iter_fasta(panel_fasta), start=1):
        original_id = parse_header_id(header)
        new_id = f"SHUYU_{idx:06d}_{safe_id(original_id)}"
        category = infer_panel_category(header)
        output_handle.write(f">{new_id} {header[1:].strip()}\n")
        for chunk in chunks:
            if chunk:
                output_handle.write(chunk.upper() + "\n")
        row = {
            "reference_id": new_id,
            "category": category,
            "source": f"Shuyu_masked_panel:{panel_fasta}",
            "description": header[1:].strip(),
            "original_reference_id": original_id,
            "original_header": header[1:].strip(),
        }
        map_rows.append(row)
        inventory_rows.append(row)


def write_table(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Build hg38/HERV/LINE1 plus Shuyu masked viral panel reference.")
    parser.add_argument("--base-fasta", type=Path, required=True)
    parser.add_argument("--base-map", type=Path, required=True)
    parser.add_argument("--panel-fasta", type=Path, required=True)
    parser.add_argument("--output-fasta", type=Path, required=True)
    parser.add_argument("--output-map", type=Path, required=True)
    parser.add_argument("--panel-inventory", type=Path, required=True)
    parser.add_argument("--base-category", action="append", default=DEFAULT_BASE_CATEGORIES)
    parser.add_argument("--index", action="store_true")
    args = parser.parse_args()

    for path in (args.base_fasta, args.base_map, args.panel_fasta):
        if not path.exists():
            raise FileNotFoundError(path)

    base_map = load_reference_map(args.base_map)
    base_categories = set(args.base_category)
    map_rows: list[dict[str, str]] = []
    inventory_rows: list[dict[str, str]] = []
    columns = ["reference_id", "category", "source", "description", "original_reference_id", "original_header"]

    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with args.output_fasta.open("w", encoding="utf-8") as out:
        append_base_records(args.base_fasta, base_map, base_categories, out, map_rows)
        append_panel_records(args.panel_fasta, out, map_rows, inventory_rows)

    write_table(args.output_map, map_rows, columns)
    write_table(args.panel_inventory, inventory_rows, columns)

    if args.index:
        subprocess.run(["bwa", "index", str(args.output_fasta)], check=True)

    print(f"Wrote {args.output_fasta}")
    print(f"Wrote {args.output_map}")
    print(f"Wrote {args.panel_inventory}")
    print(f"Reference records: {len(map_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
