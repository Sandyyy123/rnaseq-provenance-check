"""Provenance and consistency audit for a finished bulk RNA-seq run.

Reads a run manifest (JSON), the sample sheet and the gene-count matrix, and reports
findings in five areas: sample identity, file integrity, reference consistency,
software pinning and QC thresholds. Standard library only, so it runs anywhere.
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
import statistics
from dataclasses import dataclass, asdict
from pathlib import Path

FLOATING_VERSIONS = {"", "latest", "stable", "main", "master", "dev", "*"}


@dataclass
class Finding:
    area: str
    severity: str  # "error" | "warning" | "ok"
    message: str


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 16), b""):
            h.update(block)
    return h.hexdigest()


def read_sample_sheet(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def read_counts(path: Path) -> tuple[list[str], dict[str, list[int]]]:
    """Return (sample columns, {gene_id: counts})."""
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh))
    header, body = rows[0], rows[1:]
    samples = header[1:]
    return samples, {r[0]: [int(x) for x in r[1:]] for r in body}


def check_samples(sheet: list[dict], fastq_entries: list[dict], count_samples: list[str]) -> list[Finding]:
    out = []
    ids = [r["sample_id"] for r in sheet]
    dupes = sorted({s for s in ids if ids.count(s) > 1})
    if dupes:
        out.append(Finding("samples", "error", f"duplicate sample ids in sample sheet: {', '.join(dupes)}"))
    sheet_set, fq_set, cnt_set = set(ids), {f["sample_id"] for f in fastq_entries}, set(count_samples)
    for name, other in (("FASTQ list", fq_set), ("count matrix", cnt_set)):
        missing, extra = sorted(sheet_set - other), sorted(other - sheet_set)
        if missing:
            out.append(Finding("samples", "error", f"in sample sheet but not in {name}: {', '.join(missing)}"))
        if extra:
            out.append(Finding("samples", "error", f"in {name} but not in sample sheet: {', '.join(extra)}"))
    if not out:
        out.append(Finding("samples", "ok", f"{len(ids)} samples consistent across sample sheet, FASTQ list and count matrix"))
    return out


def check_files(base: Path, fastq_entries: list[dict]) -> list[Finding]:
    out = []
    for f in fastq_entries:
        p = base / f["path"]
        if not p.is_file():
            out.append(Finding("files", "error", f"missing file {f['path']}"))
        elif sha256(p) != f["sha256"]:
            out.append(Finding("files", "error", f"checksum mismatch for {f['path']} (file changed after the manifest was written)"))
    if not out:
        out.append(Finding("files", "ok", f"{len(fastq_entries)} FASTQ checksums match the manifest"))
    return out


def check_reference(manifest: dict, gene_ids: list[str]) -> list[Finding]:
    out = []
    steps = manifest["steps"]
    builds = {s["name"]: s.get("genome_build") for s in steps if s.get("genome_build")}
    annots = {s["name"]: s.get("annotation") for s in steps if s.get("annotation")}
    if len(set(builds.values())) > 1:
        out.append(Finding("reference", "error", "genome build differs between steps: " + "; ".join(f"{k}={v}" for k, v in builds.items())))
    if len(set(annots.values())) > 1:
        out.append(Finding("reference", "error", "annotation release differs between steps: " + "; ".join(f"{k}={v}" for k, v in annots.items())))
    versioned = sum(1 for g in gene_ids if re.fullmatch(r"ENSG\d+\.\d+", g))
    if 0 < versioned < len(gene_ids):
        out.append(Finding("reference", "warning", f"gene ids mix versioned and unversioned Ensembl ids ({versioned} of {len(gene_ids)} versioned)"))
    if not out:
        out.append(Finding("reference", "ok", "one genome build and one annotation release used throughout"))
    return out


def check_software(manifest: dict) -> list[Finding]:
    out = []
    for s in manifest["steps"]:
        v = str(s.get("version", "")).strip().lower()
        if v in FLOATING_VERSIONS:
            out.append(Finding("software", "error", f"{s['name']} ({s['tool']}) has no pinned version ('{s.get('version', '')}')"))
        if not s.get("container_digest") and not s.get("conda_lock"):
            out.append(Finding("software", "warning", f"{s['name']} has no container digest or lock file"))
    if not out:
        out.append(Finding("software", "ok", "every step pinned to a version and an environment"))
    return out


def check_qc(manifest: dict, count_samples: list[str], counts: dict[str, list[int]]) -> list[Finding]:
    out = []
    th = manifest["qc_thresholds"]
    for sid, m in manifest["qc_metrics"].items():
        if m["mapping_rate"] < th["min_mapping_rate"]:
            out.append(Finding("qc", "warning", f"{sid}: mapping rate {m['mapping_rate']:.0%} below {th['min_mapping_rate']:.0%}"))
        if m["rrna_fraction"] > th["max_rrna_fraction"]:
            out.append(Finding("qc", "warning", f"{sid}: rRNA fraction {m['rrna_fraction']:.0%} above {th['max_rrna_fraction']:.0%}"))
    libsizes = {s: sum(c[i] for c in counts.values()) for i, s in enumerate(count_samples)}
    med = statistics.median(libsizes.values())
    for s, n in libsizes.items():
        if n < th["min_library_fraction_of_median"] * med:
            out.append(Finding("qc", "warning", f"{s}: assigned counts {n:,} are under {th['min_library_fraction_of_median']:.0%} of the median ({med:,.0f})"))
    if not out:
        out.append(Finding("qc", "ok", "all samples within QC thresholds"))
    return out


def audit(run_dir: Path) -> dict:
    run_dir = Path(run_dir)
    manifest = json.loads((run_dir / "manifest.json").read_text())
    sheet = read_sample_sheet(run_dir / manifest["sample_sheet"])
    count_samples, counts = read_counts(run_dir / manifest["counts"])
    findings = (check_samples(sheet, manifest["fastq"], count_samples)
                + check_files(run_dir, manifest["fastq"])
                + check_reference(manifest, list(counts))
                + check_software(manifest)
                + check_qc(manifest, count_samples, counts))
    inputs = sorted([manifest["sample_sheet"], manifest["counts"], "manifest.json"] + [f["path"] for f in manifest["fastq"] if (run_dir / f["path"]).is_file()])
    return {
        "run": manifest["run_id"],
        "n_samples": len(sheet),
        "n_genes": len(counts),
        "summary": {sev: sum(1 for f in findings if f.severity == sev) for sev in ("error", "warning", "ok")},
        "findings": [asdict(f) for f in findings],
        "input_hashes": {p: sha256(run_dir / p) for p in inputs},
    }
