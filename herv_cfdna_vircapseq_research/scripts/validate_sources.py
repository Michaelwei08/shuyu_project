from __future__ import annotations

import argparse
import csv
import json
import re
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable
from urllib.parse import urlparse

import requests
import yaml
from bs4 import BeautifulSoup
from requests import Response, Session


URL_RE = re.compile(r"https?://[^\s)\]`\"'>]+")
TEXT_FILE_SUFFIXES = {".md", ".json", ".yaml", ".yml", ".txt"}
USER_AGENT = "herv-cfdna-research-validator/1.0"


@dataclass
class SourceRecord:
    source_url: str
    normalized_url: str
    domain: str
    category: str
    status_code: int | None = None
    ok: bool = False
    final_url: str | None = None
    content_type: str | None = None
    title: str | None = None
    warning: str | None = None
    repo_full_name: str | None = None
    repo_stars: int | None = None
    repo_updated_at: str | None = None
    pubmed_id: str | None = None
    pmcid: str | None = None
    doi: str | None = None
    publication_year: str | None = None


def normalize_url(url: str) -> str:
    url = url.rstrip(".,;")
    parsed = urlparse(url)
    scheme = parsed.scheme.lower() or "https"
    netloc = parsed.netloc.lower()
    if netloc == "www.dfam.org":
        netloc = "dfam.org"
    if netloc == "www.repeatmasker.org":
        netloc = "repeatmasker.org"
    path = parsed.path.rstrip("/") or "/"
    if parsed.query:
        return f"{scheme}://{netloc}{path}?{parsed.query}"
    return f"{scheme}://{netloc}{path}"


def classify_url(url: str) -> str:
    parsed = urlparse(url)
    domain = parsed.netloc.lower()
    if domain.endswith("doi.org"):
        return "doi"
    if "github.com" in domain:
        return "github"
    if "pubmed.ncbi.nlm.nih.gov" in domain:
        return "pubmed"
    if "pmc.ncbi.nlm.nih.gov" in domain:
        return "pmc"
    if path_looks_like_pdf(parsed.path):
        return "pdf"
    return "web"


def path_looks_like_pdf(path: str) -> bool:
    return path.lower().endswith(".pdf")


def is_within(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


def iter_project_urls(root: Path, excluded_dirs: Iterable[Path]) -> Iterable[str]:
    seen: set[str] = set()
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        if any(is_within(path, excluded_dir) for excluded_dir in excluded_dirs):
            continue
        if path.suffix.lower() not in TEXT_FILE_SUFFIXES:
            continue
        text = path.read_text(encoding="utf-8")
        for match in URL_RE.findall(text):
            normalized = normalize_url(match)
            if normalized not in seen:
                seen.add(normalized)
                yield normalized


def request_with_fallback(session: Session, url: str) -> Response:
    response = session.get(
        url,
        timeout=30,
        allow_redirects=True,
        headers={"User-Agent": USER_AGENT},
    )
    return response


def parse_title_from_response(response: Response) -> str | None:
    content_type = response.headers.get("content-type", "")
    if "html" not in content_type.lower():
        return None
    soup = BeautifulSoup(response.text, "html.parser")
    for attr in ("og:title", "twitter:title"):
        tag = soup.find("meta", attrs={"property": attr}) or soup.find(
            "meta", attrs={"name": attr}
        )
        if tag and tag.get("content"):
            return " ".join(tag["content"].split())
    if soup.title and soup.title.string:
        return " ".join(soup.title.string.split())
    return None


def looks_like_placeholder_title(title: str | None) -> bool:
    if not title:
        return True
    normalized = title.lower()
    return (
        "recaptcha" in normalized
        or "checking your browser" in normalized
        or "redirecting" in normalized
        or "just a moment" in normalized
    )


def github_metadata(session: Session, url: str) -> dict[str, str | int | None]:
    parsed = urlparse(url)
    parts = [part for part in parsed.path.split("/") if part]
    if len(parts) < 2:
        return {}
    owner, repo = parts[0], parts[1]
    api_url = f"https://api.github.com/repos/{owner}/{repo}"
    try:
        response = session.get(api_url, timeout=30, headers={"User-Agent": USER_AGENT})
    except requests.RequestException:
        return {}
    if not response.ok:
        return {}
    payload = response.json()
    return {
        "repo_full_name": payload.get("full_name"),
        "repo_stars": payload.get("stargazers_count"),
        "repo_updated_at": payload.get("updated_at"),
        "title": payload.get("description"),
    }


def crossref_metadata(session: Session, doi: str | None) -> dict[str, str | None]:
    if not doi:
        return {}
    api_url = f"https://api.crossref.org/works/{doi}"
    try:
        response = session.get(api_url, timeout=30, headers={"User-Agent": USER_AGENT})
    except requests.RequestException:
        return {}
    if not response.ok:
        return {}
    message = response.json().get("message", {})
    title_list = message.get("title") or []
    year = None
    issued = message.get("issued", {}).get("date-parts", [])
    if issued and issued[0]:
        year = str(issued[0][0])
    return {
        "title": title_list[0] if title_list else None,
        "publication_year": year,
        "doi": message.get("DOI"),
    }


def parse_pubmed_id(url: str) -> str | None:
    match = re.search(r"/(\d+)/?$", urlparse(url).path)
    return match.group(1) if match else None


def parse_pmcid(url: str) -> str | None:
    match = re.search(r"/articles/(PMC\d+)/?$", urlparse(url).path)
    return match.group(1) if match else None


def pubmed_metadata(session: Session, pmid: str | None) -> dict[str, str | None]:
    if not pmid:
        return {}
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=pubmed&id={pmid}&retmode=xml"
    )
    try:
        response = session.get(url, timeout=30, headers={"User-Agent": USER_AGENT})
    except requests.RequestException:
        return {}
    if not response.ok:
        return {}
    soup = BeautifulSoup(response.text, "xml")
    title_tag = soup.find("ArticleTitle")
    year_tag = soup.find("PubDate").find("Year") if soup.find("PubDate") else None
    article_id_tags = soup.find_all("ArticleId")
    doi = None
    for tag in article_id_tags:
        if tag.get("IdType") == "doi":
            doi = tag.text.strip()
            break
    return {
        "title": title_tag.get_text(" ", strip=True) if title_tag else None,
        "publication_year": year_tag.text.strip() if year_tag else None,
        "doi": doi,
    }


def pmc_to_pubmed(session: Session, pmcid: str | None) -> dict[str, str | None]:
    if not pmcid:
        return {}
    url = (
        "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
        f"?ids={pmcid}&format=json"
    )
    try:
        response = session.get(url, timeout=30, headers={"User-Agent": USER_AGENT})
    except requests.RequestException:
        return {}
    if not response.ok:
        return {}
    payload = response.json()
    records = payload.get("records") or []
    if not records:
        return {}
    record = records[0]
    return {
        "pmcid": record.get("pmcid"),
        "pubmed_id": record.get("pmid"),
        "doi": record.get("doi"),
    }


def validate_url(session: Session, url: str) -> SourceRecord:
    category = classify_url(url)
    parsed = urlparse(url)
    record = SourceRecord(
        source_url=url,
        normalized_url=normalize_url(url),
        domain=parsed.netloc.lower(),
        category=category,
    )

    try:
        response = request_with_fallback(session, url)
        record.status_code = response.status_code
        record.ok = response.ok
        record.final_url = response.url
        record.content_type = response.headers.get("content-type")
        record.title = parse_title_from_response(response)
        lower_text = response.text[:2000].lower()
        if "recaptcha" in lower_text:
            record.warning = "blocked_by_recaptcha"
        elif response.status_code >= 400:
            record.warning = f"http_{response.status_code}"
    except requests.RequestException as exc:
        record.warning = type(exc).__name__
        return record

    if category == "doi":
        doi = urlparse(url).path.lstrip("/")
        meta = crossref_metadata(session, doi)
        for key, value in meta.items():
            if looks_like_placeholder_title(record.title) and key == "title" and value is not None:
                setattr(record, key, value)
                continue
            if getattr(record, key) is None and value is not None:
                setattr(record, key, value)
    elif category == "github":
        meta = github_metadata(session, url)
        for key, value in meta.items():
            setattr(record, key, value)
    elif category == "pubmed":
        record.pubmed_id = parse_pubmed_id(url)
        meta = pubmed_metadata(session, record.pubmed_id)
        for key, value in meta.items():
            if getattr(record, key) is None and value is not None:
                setattr(record, key, value)
    elif category == "pmc":
        record.pmcid = parse_pmcid(url)
        ids = pmc_to_pubmed(session, record.pmcid)
        for key, value in ids.items():
            if value is not None:
                setattr(record, key, value)
        meta = pubmed_metadata(session, record.pubmed_id)
        for key, value in meta.items():
            if key == "title":
                if looks_like_placeholder_title(record.title) and value is not None:
                    setattr(record, key, value)
                continue
            if getattr(record, key) is None and value is not None:
                setattr(record, key, value)
        if record.warning is None and record.title and "recaptcha" in record.title.lower():
            record.warning = "blocked_by_recaptcha"
    elif record.doi is None:
        doi_match = re.search(r"/(10\.\d{4,9}/[^?#]+)", url)
        if doi_match:
            record.doi = doi_match.group(1)
            meta = crossref_metadata(session, record.doi)
            for key, value in meta.items():
                if looks_like_placeholder_title(record.title) and key == "title" and value is not None:
                    setattr(record, key, value)
                    continue
                if getattr(record, key) is None and value is not None:
                    setattr(record, key, value)

    return record


def markdown_report(records: list[SourceRecord]) -> str:
    total = len(records)
    ok_count = sum(1 for record in records if record.ok)
    blocked = [record for record in records if record.warning]
    github = [record for record in records if record.category == "github"]
    lines = [
        "# Source Validation",
        "",
        f"- Total sources: {total}",
        f"- HTTP OK: {ok_count}",
        f"- Sources with warnings: {len(blocked)}",
        f"- GitHub repositories checked: {len(github)}",
        "",
        "## Warnings",
        "",
    ]
    if blocked:
        for record in blocked:
            lines.append(
                f"- `{record.source_url}`: `{record.warning}`"
                + (f" (status {record.status_code})" if record.status_code else "")
            )
    else:
        lines.append("- None")

    lines.extend(["", "## Source Table", "", "| Domain | Title | Status | Notes |", "| --- | --- | --- | --- |"])
    for record in sorted(records, key=lambda item: (item.domain, item.source_url)):
        title = record.title or record.repo_full_name or ""
        status = str(record.status_code) if record.status_code is not None else "n/a"
        notes: list[str] = []
        if record.warning:
            notes.append(record.warning)
        if record.repo_updated_at:
            notes.append(f"repo updated {record.repo_updated_at[:10]}")
        if record.pubmed_id:
            notes.append(f"PMID {record.pubmed_id}")
        if record.pmcid:
            notes.append(record.pmcid)
        lines.append(
            f"| {record.domain} | {title.replace('|', '/')} | {status} | {'; '.join(notes)} |"
        )
    return "\n".join(lines) + "\n"


def write_outputs(records: list[SourceRecord], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "source_validation.json"
    csv_path = output_dir / "source_validation.csv"
    md_path = output_dir / "source_validation.md"
    summary_path = output_dir / "source_summary.yaml"

    json_path.write_text(
        json.dumps([asdict(record) for record in records], indent=2),
        encoding="utf-8",
    )

    fieldnames = list(asdict(records[0]).keys()) if records else []
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for record in records:
            writer.writerow(asdict(record))

    md_path.write_text(markdown_report(records), encoding="utf-8")

    summary = {
        "generated_at_unix": int(time.time()),
        "total_sources": len(records),
        "ok_sources": sum(1 for record in records if record.ok),
        "sources_with_warnings": sum(1 for record in records if record.warning),
        "domains": sorted({record.domain for record in records}),
    }
    summary_path.write_text(yaml.safe_dump(summary, sort_keys=False), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-root",
        default=".",
        help="Project directory containing markdown/json/yaml research artifacts.",
    )
    parser.add_argument(
        "--output-dir",
        default="validation",
        help="Directory where validation outputs will be written.",
    )
    args = parser.parse_args()

    project_root = Path(args.project_root).resolve()
    output_dir = Path(args.output_dir).resolve()

    excluded_dirs = [output_dir]
    urls = sorted(iter_project_urls(project_root, excluded_dirs))
    session = requests.Session()
    records = [validate_url(session, url) for url in urls]
    write_outputs(records, output_dir)


if __name__ == "__main__":
    main()
