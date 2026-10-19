"""Create the synthetic demo run in examples/demo_run.

The data are invented: 6 samples, 40 genes, tiny placeholder FASTQ files. Four problems are
planted on purpose so the audit has something real to find:
  1. the quantification step uses a different annotation release than the alignment index;
  2. the QC report tool is unpinned ('latest');
  3. sample S6 is in the count matrix but missing from the sample sheet;
  4. sample S4 has a low mapping rate and a small library.
One FASTQ is edited after the manifest was written, so its checksum no longer matches.
"""
import csv
import hashlib
import json
import random
from pathlib import Path

ROOT = Path(__file__).parent / "examples" / "demo_run"


def main() -> None:
    rnd = random.Random(7)
    (ROOT / "fastq").mkdir(parents=True, exist_ok=True)
    samples = ["S1", "S2", "S3", "S4", "S5"]
    conditions = ["control", "control", "control", "treated", "treated"]
    with open(ROOT / "sample_sheet.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample_id", "condition", "batch"])
        for s, c in zip(samples, conditions):
            w.writerow([s, c, "B1" if s in ("S1", "S2", "S4") else "B2"])

    fastq = []
    for s in samples:
        p = ROOT / "fastq" / f"{s}_R1.fastq"
        p.write_text(f"@{s}_read1\nACGTACGTAC\n+\nIIIIIIIIII\n")
        fastq.append({"sample_id": s, "path": f"fastq/{s}_R1.fastq", "sha256": hashlib.sha256(p.read_bytes()).hexdigest()})
    # planted: S2 file edited after the manifest was written
    (ROOT / "fastq" / "S2_R1.fastq").write_text("@S2_read1\nACGTACGTAA\n+\nIIIIIIIIII\n")

    count_samples = samples + ["S6"]
    genes = [f"ENSG{100000 + i:011d}.{rnd.randint(1, 9)}" for i in range(36)] + [f"ENSG{200000 + i:011d}" for i in range(4)]
    with open(ROOT / "counts.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["gene_id"] + count_samples)
        for g in genes:
            base = rnd.randint(50, 5000)
            row = []
            for s in count_samples:
                scale = 0.25 if s == "S4" else 1.0
                row.append(int(base * scale * rnd.uniform(0.8, 1.2)))
            w.writerow([g] + row)

    manifest = {
        "run_id": "demo-bulk-rnaseq-001",
        "sample_sheet": "sample_sheet.csv",
        "counts": "counts.csv",
        "fastq": fastq,
        "steps": [
            {"name": "align", "tool": "STAR", "version": "2.7.11b", "genome_build": "GRCh38", "annotation": "GENCODE v44",
             "container_digest": "sha256:1111aaaa"},
            {"name": "quantify", "tool": "featureCounts", "version": "2.0.6", "genome_build": "GRCh38", "annotation": "GENCODE v43",
             "container_digest": "sha256:2222bbbb"},
            {"name": "qc_report", "tool": "MultiQC", "version": "latest"},
        ],
        "qc_thresholds": {"min_mapping_rate": 0.70, "max_rrna_fraction": 0.10, "min_library_fraction_of_median": 0.5},
        "qc_metrics": {
            "S1": {"mapping_rate": 0.91, "rrna_fraction": 0.03},
            "S2": {"mapping_rate": 0.89, "rrna_fraction": 0.04},
            "S3": {"mapping_rate": 0.90, "rrna_fraction": 0.02},
            "S4": {"mapping_rate": 0.58, "rrna_fraction": 0.14},
            "S5": {"mapping_rate": 0.88, "rrna_fraction": 0.05},
        },
    }
    (ROOT / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"demo run written to {ROOT}")


if __name__ == "__main__":
    main()
