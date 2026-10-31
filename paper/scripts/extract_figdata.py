"""Compact figure inputs from the banhbo benchmark outputs (run on banhbo, results copied to the Mac).
timing_<set>.tsv: per-locus seconds for every tool; pip_<set>.tsv: PIP pairs vs susieR 0.14.2 for variants with
PIP >= 0.001 in either tool (the rest are near zero and only crowd the plot), plus per-locus max |difference|;
coloc_pairs.tsv: tested pairs (PP.H4 not missing) for flashLocus and R; ld_times.tsv copied as is."""
import os, sys, numpy as np, pandas as pd
B = os.environ.get("BENCH_ROOT", "."); O = sys.argv[1]; os.makedirs(O, exist_ok=True)
tools = ["susier014", "susier016", "pysusie", "pysusie_exact", "rs_sym", "rs_auto"]
for s, man in [("sim", f"{B}/bench/sim2/manifest.tsv"), ("real", f"{B}/real/all.tsv")]:
    ids = [l.split("\t")[0] for l in open(man) if l.strip()]
    T = pd.concat({t: pd.read_csv(f"{B}/final/{s}/{t}/timing.tsv", sep="\t").set_index("id")[["p", "seconds"]] for t in tools}, axis=1)
    t = pd.DataFrame({"id": T.index, "p": T[("susier014", "p")].values})
    for k in tools: t[k] = T[(k, "seconds")].values
    t.to_csv(f"{O}/timing_{s}.tsv", sep="\t", index=False)
    rows, mx = [], []
    for i in ids:
        ref = np.loadtxt(f"{B}/final/{s}/susier014/{i}.pip")
        d = {"id": i}
        for k in ["rs_sym", "rs_auto", "pysusie", "susier016"]:
            x = np.loadtxt(f"{B}/final/{s}/{k}/{i}.pip"); d[k] = np.abs(x - ref).max()
            keep = (x >= 1e-3) | (ref >= 1e-3)
            rows += [(i, k, a, b) for a, b in zip(ref[keep], x[keep])]
        mx.append(d)
    pd.DataFrame(rows, columns=["id", "tool", "pip_susier014", "pip_tool"]).to_csv(f"{O}/pip_{s}.tsv.gz", sep="\t", index=False, float_format="%.6g")
    pd.DataFrame(mx).to_csv(f"{O}/pipmaxdiff_{s}.tsv", sep="\t", index=False, float_format="%.3e")
P = []
for n in ["uterus", "ovary", "vagina", "adipose_visceral", "blood", "adipose_subcut"]:
    a = pd.read_csv(f"{B}/atlas/full/{n}.coloc.tsv", sep="\t"); b = pd.read_csv(f"{B}/atlas/full/{n}.R.tsv", sep="\t")
    m = a.merge(b, on=["locus", "tissue", "gene", "gwas_L", "qtl_L"], suffixes=("_fl", "_r"))
    m = m[m["PP.H4_r"].notna()]
    P.append(m[["locus", "tissue", "gene", "gwas_L", "qtl_L", "PP.H4_fl", "PP.H4_r"]])
pd.concat(P).to_csv(f"{O}/coloc_pairs.tsv.gz", sep="\t", index=False, float_format="%.17g")
os.system(f"cp {B}/ldbench/ld_times.tsv {B}/atlas/full_v2seq/timing.tsv {O}/ 2>/dev/null; cp {B}/atlas/full_v2seq/timing.tsv {O}/atlas_flashlocus_timing.tsv")
os.system(f"grep -h -E 'R atlas|wall' {B}/atlas/logs/r.log {B}/atlas/logs/r_big.log > {O}/atlas_r_logs.txt; grep -h wall {B}/final/logs/*_tp.log > {O}/throughput.txt")
for s in ["sim", "real"]:
    os.system(f"for t in rs_sym rs_auto; do echo \"$t $(cat {B}/final/{s}/tp_$t/batch.log 2>/dev/null | tail -1)\"; done >> {O}/throughput_{s}_rs.txt")
print("done")
