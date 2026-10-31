"""Python API regression tests against susieR 0.14.2 and coloc 5.2.3 fixtures (tests/data)."""
import csv, math, os, subprocess, sys
import numpy as np
import flashlocus

D = os.path.join(os.path.dirname(__file__), "..", "data")


def test_susie_rss_matches_susier():
    z = np.loadtxt(os.path.join(D, "sim000.z"))
    R = np.fromfile(os.path.join(D, "sim000.ld.bin")).reshape(len(z), len(z))
    fit = flashlocus.susie_rss(z, R, 100000)
    ref = np.loadtxt(os.path.join(D, "sim000.susier014.pip"))
    assert np.abs(fit["pip"] - ref).max() < 1e-6
    cs = [l.split("\t") for l in open(os.path.join(D, "sim000.susier014.cs")) if l.strip()]
    assert len(fit["sets"]) == len(cs)
    assert [m + 1 for m in fit["sets"][0]["members"]] == [int(x) for x in cs[0][1].split(",")]
    assert fit["alpha"].shape == (10, len(z)) and fit["converged"]


def _bf(path):
    rows = list(csv.reader(open(path), delimiter="\t"))[1:]
    return np.array([[float(x) for x in r[1:]] for r in rows]).T, [r[0] for r in rows]


def test_coloc_matches_r():
    for case in ["partial", "dropall"]:
        b1, s1 = _bf(os.path.join(D, f"coloc_{case}_bf1.tsv"))
        b2, s2 = _bf(os.path.join(D, f"coloc_{case}_bf2.tsv"))
        got = {(r["idx1"], r["idx2"]): r for r in flashlocus.coloc_bf_bf(b1, s1, b2, s2)}
        exp = list(csv.DictReader(open(os.path.join(D, f"coloc_{case}_R.tsv")), delimiter="\t"))
        assert len(got) == len(exp)
        for e in exp:
            g = got[(int(e["idx1"]) - 1, int(e["idx2"]) - 1)]
            assert g["hit1"] == e["hit1"] and g["hit2"] == e["hit2"] and g["nsnps"] == int(e["nsnps"])
            for h in range(5):
                v = float(e[f"PP.H{h}.abf"])
                assert (math.isnan(v) and math.isnan(g[f"PP.H{h}"])) or abs(g[f"PP.H{h}"] - v) < 1e-10


def test_cli_entry_point():
    r = subprocess.run([sys.executable, "-m", "flashlocus"], capture_output=True, text=True)
    assert "usage: flashlocus" in r.stderr
