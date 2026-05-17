from __future__ import annotations

import argparse
import csv
import subprocess
import tempfile
from pathlib import Path


DEFAULT_NCBI_RECORDS = [
    ("K03455", "HIV1", "Human immunodeficiency virus 1 HXB2"),
    ("M15390", "HIV2", "Human immunodeficiency virus 2 ROD"),
    ("J02029", "HTLV1", "Human T-cell lymphotropic virus type 1 ATK"),
    ("M10060", "HTLV2", "Human T-cell lymphotropic virus type 2 Mo"),
    ("AY037928", "HERV", "Human endogenous retrovirus K113"),
    ("L19088", "LINE1", "Human LINE-1 L1.3"),
]


def parse_header_id(header: str) -> str:
    return header[1:].strip().split()[0]


def iter_fasta_records(path: Path):
    header = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, chunks
                header = line
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        yield header, chunks


def append_fasta(
    source_fasta: Path,
    category: str,
    output_handle,
    map_rows: list[dict[str, str]],
    source_label: str,
) -> None:
    for header, chunks in iter_fasta_records(source_fasta):
        output_handle.write(header + "\n")
        for chunk in chunks:
            if chunk:
                output_handle.write(chunk + "\n")
        map_rows.append(
            {
                "reference_id": parse_header_id(header),
                "category": category,
                "source": source_label,
                "description": header[1:].strip(),
            }
        )


def efetch_accession(accession: str, output_path: Path) -> None:
    with output_path.open("w", encoding="utf-8") as handle:
        subprocess.run(
            ["efetch", "-db", "nucleotide", "-format", "fasta", "-id", accession],
            stdout=handle,
            check=True,
        )


def parse_extra_fasta(value: str) -> tuple[str, Path]:
    if ":" not in value:
        raise ValueError("--extra-fasta must be CATEGORY:/path/to/file.fa")
    category, path_text = value.split(":", 1)
    category = category.strip()
    if not category:
        raise ValueError("extra FASTA category cannot be empty")
    path = Path(path_text)
    if not path.exists():
        raise FileNotFoundError(path)
    return category, path


def write_map(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["reference_id", "category", "source", "description"])
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build a competitive retrovirus/HERV/LINE1 FASTA and reference-category map."
    )
    parser.add_argument("--output-fasta", type=Path, required=True)
    parser.add_argument("--output-map", type=Path, required=True)
    parser.add_argument(
        "--from-efetch-defaults",
        action="store_true",
        help="Fetch HIV1, HIV2, HTLV1, HTLV2, HERV-K113, and LINE1 L1.3 from NCBI using efetch.",
    )
    parser.add_argument(
        "--extra-fasta",
        action="append",
        default=[],
        help="Append records from a local FASTA as CATEGORY:/path/to/file.fa. Can be repeated.",
    )
    parser.add_argument("--index", action="store_true", help="Run bwa index on the output FASTA.")
    args = parser.parse_args()

    if not args.from_efetch_defaults and not args.extra_fasta:
        raise SystemExit("Provide --from-efetch-defaults and/or one or more --extra-fasta CATEGORY:path entries.")

    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    map_rows: list[dict[str, str]] = []

    with args.output_fasta.open("w", encoding="utf-8") as out:
        if args.from_efetch_defaults:
            with tempfile.TemporaryDirectory() as tmpdir:
                tmp = Path(tmpdir)
                for accession, category, description in DEFAULT_NCBI_RECORDS:
                    fetched = tmp / f"{accession}.fa"
                    efetch_accession(accession, fetched)
                    append_fasta(fetched, category, out, map_rows, f"NCBI:{accession}:{description}")

        for item in args.extra_fasta:
            category, fasta_path = parse_extra_fasta(item)
            append_fasta(fasta_path, category, out, map_rows, f"local:{fasta_path}")

    write_map(args.output_map, map_rows)

    if args.index:
        subprocess.run(["bwa", "index", str(args.output_fasta)], check=True)

    print(f"Wrote {args.output_fasta}")
    print(f"Wrote {args.output_map}")
    print(f"Reference records: {len(map_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
