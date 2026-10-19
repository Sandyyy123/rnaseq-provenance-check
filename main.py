"""Run the provenance check on a run directory and print the findings.

Usage: python main.py [RUN_DIR] [--json OUT.json]
Default RUN_DIR is the synthetic demo in examples/demo_run.
"""
import argparse
import json
from pathlib import Path

from rnaprov.audit import audit

ICON = {"error": "ERROR  ", "warning": "WARNING", "ok": "OK     "}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir", nargs="?", default=str(Path(__file__).parent / "examples" / "demo_run"))
    ap.add_argument("--json", dest="json_out")
    args = ap.parse_args()
    report = audit(Path(args.run_dir))
    print(f"run {report['run']}: {report['n_samples']} samples in sheet, {report['n_genes']} genes")
    for f in report["findings"]:
        print(f"{ICON[f['severity']]} [{f['area']}] {f['message']}")
    print("summary:", report["summary"])
    if args.json_out:
        Path(args.json_out).write_text(json.dumps(report, indent=2))
    return 1 if report["summary"]["error"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
