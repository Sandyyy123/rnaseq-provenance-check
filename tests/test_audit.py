import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from rnaprov.audit import audit, check_reference, check_software  # noqa: E402


@pytest.fixture()
def demo(tmp_path):
    subprocess.run([sys.executable, str(ROOT / "make_demo.py")], check=True, capture_output=True)
    dst = tmp_path / "run"
    shutil.copytree(ROOT / "examples" / "demo_run", dst)
    return dst


def messages(report, area, severity):
    return [f["message"] for f in report["findings"] if f["area"] == area and f["severity"] == severity]


def test_demo_finds_all_planted_problems(demo):
    r = audit(demo)
    assert r["summary"]["error"] == 4
    assert any("S6" in m for m in messages(r, "samples", "error"))
    assert any("S2_R1.fastq" in m for m in messages(r, "files", "error"))
    assert any("annotation release differs" in m for m in messages(r, "reference", "error"))
    assert any("MultiQC" in m for m in messages(r, "software", "error"))


def test_qc_flags_s4_only(demo):
    r = audit(demo)
    qc = messages(r, "qc", "warning")
    assert qc and all(m.startswith("S4") for m in qc)
    assert len(qc) == 3  # mapping rate, rRNA, library size


def test_mixed_gene_id_versions_warned(demo):
    r = audit(demo)
    assert any("mix versioned" in m for m in messages(r, "reference", "warning"))


def test_clean_run_passes(demo):
    m = json.loads((demo / "manifest.json").read_text())
    m["steps"][1]["annotation"] = "GENCODE v44"
    m["steps"][2]["version"] = "1.25"
    m["steps"][2]["conda_lock"] = "env.lock"
    (demo / "manifest.json").write_text(json.dumps(m))
    assert not check_reference(m, ["ENSG00000000001.1"])[0].severity == "error"
    assert all(f.severity == "ok" for f in check_software(m))


def test_genome_build_mismatch_is_error():
    m = {"steps": [{"name": "a", "genome_build": "GRCh38"}, {"name": "b", "genome_build": "GRCh37"}]}
    f = check_reference(m, ["ENSG00000000001"])
    assert f[0].severity == "error" and "genome build" in f[0].message


def test_input_hashes_recorded(demo):
    r = audit(demo)
    assert "manifest.json" in r["input_hashes"] and len(r["input_hashes"]["counts.csv"]) == 64


def test_cli_exit_code_signals_errors(demo):
    p = subprocess.run([sys.executable, str(ROOT / "main.py"), str(demo)], capture_output=True, text=True)
    assert p.returncode == 1 and "summary:" in p.stdout
