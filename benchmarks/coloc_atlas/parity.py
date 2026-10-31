"""Gate 1: susiers atlas vs R coloc.bf_bf, pair by pair. usage: parity.py rust.tsv r.tsv"""
import sys, numpy as np, pandas as pd
k = ["locus", "tissue", "gene", "gwas_L", "qtl_L"]
a, b = (pd.read_csv(f, sep="\t") for f in sys.argv[1:3])
m = a.merge(b, on=k, how="outer", suffixes=("_rs", "_r"), indicator=True)
print("pairs rust", len(a), "R", len(b), "| only rust", (m._merge == "left_only").sum(), "only R", (m._merge == "right_only").sum())
m = m[m._merge == "both"]
print("nsnps equal", (m.nsnps_rs == m.nsnps_r).all(), "| hit1 equal", (m.hit1_rs == m.hit1_r).mean(), "| hit2 equal", (m.hit2_rs == m.hit2_r).mean())
for h in ["PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"]:
    x, y = m[h + "_rs"].values, m[h + "_r"].values
    nan_same = (np.isnan(x) == np.isnan(y)).all()
    d = np.nanmax(np.abs(x - y)) if (~np.isnan(x)).any() else 0
    print(f"{h}: NaN pattern identical {nan_same}, max |diff| {d:.2e}")
v = m[~np.isnan(m["PP.H4_r"])]
print("tested pairs", len(v), "| H4>=0.8 rust", (v["PP.H4_rs"] >= 0.8).sum(), "R", (v["PP.H4_r"] >= 0.8).sum())
