from __future__ import annotations

import argparse
import csv
from pathlib import Path


def load_metadata(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if "sample_id" not in (reader.fieldnames or []):
            raise ValueError("Metadata must contain sample_id")
        return {row["sample_id"]: row for row in reader}


def join(viral_csv: Path, metadata_csv: Path, output_csv: Path) -> None:
    metadata = load_metadata(metadata_csv)
    with viral_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        viral_fields = list(reader.fieldnames or [])
        metadata_fields = [field for field in next(iter(metadata.values()), {}).keys() if field != "sample_id"]
        handle.seek(0)
        reader = csv.DictReader(handle)
        with output_csv.open("w", newline="", encoding="utf-8") as out:
            writer = csv.DictWriter(out, fieldnames=[*viral_fields, *metadata_fields, "metadata_join_status"])
            writer.writeheader()
            for row in reader:
                meta = metadata.get(row["sample_id"])
                output = dict(row)
                if meta:
                    for field in metadata_fields:
                        output[field] = meta.get(field, "")
                    output["metadata_join_status"] = "matched"
                else:
                    for field in metadata_fields:
                        output[field] = ""
                    output["metadata_join_status"] = "missing_metadata"
                writer.writerow(output)


def main() -> int:
    parser = argparse.ArgumentParser(description="Join viral calls to clinical metadata by sample_id.")
    parser.add_argument("viral_csv", type=Path)
    parser.add_argument("metadata_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()
    join(args.viral_csv, args.metadata_csv, args.output_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
