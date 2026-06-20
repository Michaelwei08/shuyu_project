from __future__ import annotations

import re
import subprocess
from pathlib import Path

UNMAPPED_FLAG = 0x4
SECONDARY_FLAG = 0x100
SUPPLEMENTARY_FLAG = 0x800


def samtools_exclude_flags(exclude_secondary_supplementary: bool) -> str:
    flags = UNMAPPED_FLAG
    if exclude_secondary_supplementary:
        flags |= SECONDARY_FLAG | SUPPLEMENTARY_FLAG
    return str(flags)


def parse_category_min_mapq(items: list[str]) -> dict[str, int]:
    thresholds: dict[str, int] = {}
    for item in items:
        if ":" not in item:
            raise ValueError("--category-min-mapq must be CATEGORY:MAPQ")
        category, value = item.split(":", 1)
        category = category.strip()
        if not category:
            raise ValueError("--category-min-mapq category cannot be empty")
        threshold = int(value)
        if threshold < 0:
            raise ValueError("--category-min-mapq MAPQ must be >= 0")
        thresholds[category] = threshold
    return thresholds


def aligned_len(cigar: str) -> int:
    total = 0
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        if op in "M=X":
            total += int(n)
    return total


def sam_tag(fields: list[str], tag: str) -> int | None:
    prefix = f"{tag}:i:"
    for field in fields[11:]:
        if field.startswith(prefix):
            return int(field[len(prefix) :])
    return None


def strand(flag: int, mate: bool = False) -> str:
    bit = 0x20 if mate else 0x10
    return "-" if flag & bit else "+"


def mate_reference(ref: str, mate_ref: str) -> str:
    return ref if mate_ref == "=" else mate_ref


def umi_value(read_id: str, umi_regex: str | None) -> str:
    if not umi_regex:
        raise ValueError("--dedup-mode umi requires --umi-regex")
    match = re.search(umi_regex, read_id)
    if not match:
        raise ValueError(f"UMI regex did not match read_id: {read_id}")
    if "umi" in match.groupdict():
        return match.group("umi")
    if match.groups():
        return match.group(1)
    return match.group(0)


def dedup_key(
    mode: str,
    fields: list[str],
    category: str,
    alignment_length: int,
    umi_regex: str | None,
) -> tuple[object, ...] | None:
    if mode == "none":
        return None

    read_id = fields[0]
    flag = int(fields[1])
    ref = fields[2]
    pos = int(fields[3])
    cigar = fields[5]

    if mode == "read_id":
        return (read_id, category)

    if mode == "coordinate":
        return (category, ref, pos, strand(flag), cigar, alignment_length)

    if mode == "fragment":
        mate_ref = fields[6] if len(fields) > 6 else "*"
        mate_pos = int(fields[7]) if len(fields) > 7 and fields[7].isdigit() else 0
        template_len = abs(int(fields[8])) if len(fields) > 8 and fields[8].lstrip("-").isdigit() else 0
        if mate_ref not in {"*", ""} and mate_pos > 0:
            mate_ref = mate_reference(ref, mate_ref)
            ends = sorted([(ref, pos), (mate_ref, mate_pos)])
            return (
                category,
                ends[0][0],
                ends[0][1],
                ends[1][0],
                ends[1][1],
                strand(flag),
                strand(flag, mate=True),
                template_len,
            )
        return (category, ref, pos, strand(flag), cigar, alignment_length)

    if mode == "umi":
        return (category, ref, pos, strand(flag), umi_value(read_id, umi_regex))

    raise ValueError(f"Unknown dedup mode: {mode}")


def filtered_counts(
    bam_path: Path,
    reference_map: dict[str, str],
    min_mapq: int,
    category_min_mapq: dict[str, int],
    min_aligned_length: int,
    samtools_exe: str,
    dedup_mode: str,
    umi_regex: str | None,
    require_unique_best: bool,
    exclude_secondary_supplementary: bool = False,
    target_reference_ids: set[str] | None = None,
) -> tuple[dict[str, int], dict[str, int], dict[str, int]]:
    counts: dict[str, int] = {}
    record_counts: dict[str, int] = {}
    dedup_removed: dict[str, int] = {}
    seen: set[tuple[str, ...]] = set()

    cmd = [samtools_exe, "view", "-F", samtools_exclude_flags(exclude_secondary_supplementary), str(bam_path)]
    if target_reference_ids:
        cmd.extend(sorted(target_reference_ids))

    proc = subprocess.Popen(
        cmd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert proc.stdout is not None
    try:
        for line in proc.stdout:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            ref, mapq_text, cigar = fields[2], fields[4], fields[5]
            if ref not in reference_map:
                continue
            category = reference_map[ref]
            mapq_threshold = category_min_mapq.get(category, min_mapq)
            if int(mapq_text) < mapq_threshold:
                continue
            alignment_length = aligned_len(cigar)
            if alignment_length < min_aligned_length:
                continue
            if require_unique_best:
                alignment_score = sam_tag(fields, "AS")
                suboptimal_score = sam_tag(fields, "XS")
                if alignment_score is not None and suboptimal_score is not None and alignment_score <= suboptimal_score:
                    continue

            record_counts[category] = record_counts.get(category, 0) + 1
            key = dedup_key(dedup_mode, fields, category, alignment_length, umi_regex)
            if key is not None:
                key_text = tuple(str(item) for item in key)
                if key_text in seen:
                    dedup_removed[category] = dedup_removed.get(category, 0) + 1
                    continue
                seen.add(key_text)
            counts[category] = counts.get(category, 0) + 1
    finally:
        if proc.stdout:
            proc.stdout.close()

    stderr = proc.stderr.read() if proc.stderr else ""
    if proc.stderr:
        proc.stderr.close()
    return_code = proc.wait()
    if return_code != 0:
        raise subprocess.CalledProcessError(
            return_code,
            cmd,
            stderr=stderr,
        )

    return counts, record_counts, dedup_removed
