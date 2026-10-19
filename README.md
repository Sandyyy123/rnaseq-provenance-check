> **Proprietary - All Rights Reserved.** (c) 2026 Sandeep Grover. This repository is licensed to Sandeep Grover and may **not** be used, run, copied, modified, distributed, or used to train models without prior written permission. Public visibility does not grant a license. See [LICENSE](LICENSE).

---

# rnaseq-provenance-audit

A small, dependency-free audit for a finished bulk RNA-seq run. It answers the questions a reviewer asks before trusting a count matrix: are the samples the same everywhere, are the raw files the ones that were processed, was one reference used throughout, is every tool pinned, and do any samples fail QC?

## What it checks

| Area | Check |
|---|---|
| samples | sample sheet, FASTQ list and count-matrix columns agree; no duplicate ids |
| files | each FASTQ matches the SHA-256 recorded in the run manifest |
| reference | one genome build and one annotation release across steps; Ensembl gene ids not mixed versioned/unversioned |
| software | every step has a pinned version (no `latest`) and a container digest or lock file |
| qc | mapping rate, rRNA fraction and library size against thresholds stored in the manifest |

Every report also records the SHA-256 of each input, so a later reader can confirm which files the verdict refers to.

## Run it

```bash
python make_demo.py          # writes the synthetic demo run to examples/demo_run
python main.py               # audits it; exit code 1 when errors are found
python main.py RUN_DIR --json report.json
python -m pytest -q          # 7 tests
```

The demo data are invented (5 samples in the sheet, 40 genes) with problems planted on purpose: a count column with no sample-sheet entry, a FASTQ edited after the manifest was written, alignment and quantification on different GENCODE releases, an unpinned QC tool, and one low-quality sample.

```
ERROR   [samples] in count matrix but not in sample sheet: S6
ERROR   [files] checksum mismatch for fastq/S2_R1.fastq (file changed after the manifest was written)
ERROR   [reference] annotation release differs between steps: align=GENCODE v44; quantify=GENCODE v43
ERROR   [software] qc_report (MultiQC) has no pinned version ('latest')
WARNING [qc] S4: mapping rate 58% below 70%
```

## Where it fits

In a Nextflow or Snakemake pipeline the same manifest can be generated from the workflow's own trace and parameter files, and the audit run as a final gate before results are released.
