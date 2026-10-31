"""Phase 1 gate (PLAN.md): parity vs susieR (PIP r >= 0.999 over all variants; identical CS membership in >= 95%
of loci) and single-thread speed (median per-locus time ratio susieR/Rust >= 10, pysusie/Rust >= 3)."""
import json, os, sys
import numpy as np, pandas as pd

B = os.path.join(os.environ.get("BENCH_ROOT", "."), "bench/res")
tools = ["rust", "susier", "pysusie_default", "pysusie_exact"]
T = {t: pd.read_csv(f"{B}/{t}/timing.tsv", sep="\t").set_index("id") for t in tools if os.path.exists(f"{B}/{t}/timing.tsv")}
ids = T["rust"].index


def cs_sets(t, i):
    f = f"{B}/{t}/{i}.cs.tsv" if t == "rust" else f"{B}/{t}/{i}.cs"
    if not os.path.exists(f) or os.path.getsize(f) == 0:
        return set()
    if t == "rust":
        d = pd.read_csv(f, sep="\t")
        return {tuple(sorted(int(x) for x in m.split(","))) for m in d.members.astype(str)} if len(d) else set()
    return {tuple(sorted(int(x) for x in l.split("\t")[1].split(","))) for l in open(f) if l.strip()}


res = {}
for other in [t for t in tools if t != "rust" and t in T]:
    pa, pb, same, n = [], [], 0, 0
    per = []
    for i in ids:
        fa, fb = f"{B}/rust/{i}.pip", f"{B}/{other}/{i}.pip"
        if not os.path.exists(fb):
            continue
        a, b = np.loadtxt(fa), np.loadtxt(fb)
        pa.append(a); pb.append(b)
        per.append(np.abs(a - b).max())
        s1, s2 = cs_sets("rust", i), cs_sets(other, i)
        same += int(s1 == s2); n += 1
    a, b = np.concatenate(pa), np.concatenate(pb)
    tt = T[other].reindex(ids).seconds / T["rust"].seconds
    res[other] = {"loci": n, "pip_r": float(np.corrcoef(a, b)[0, 1]), "max_abs_pip_diff": float(np.max(per)),
                  "loci_pip_diff_gt_0.01": int(np.sum(np.array(per) > 0.01)),
                  "cs_identical_share": same / n, "median_time_ratio_vs_rust": float(tt.median()),
                  "total_time_ratio": float(T[other].seconds.sum() / T["rust"].seconds.sum()),
                  "errors": int(T[other].seconds.isna().sum())}
truth = pd.read_csv(os.path.join(os.environ.get("BENCH_ROOT", "."), "bench/sim/truth.tsv"), sep="\t").set_index("id")
cov, ncs = 0, 0
for i in ids:
    c = {int(x) for x in str(truth.loc[i, "causal"]).split(",")}
    for s in cs_sets("rust", i):
        ncs += 1; cov += int(bool(set(s) & c))
res["rust_cs_containing_causal"] = cov / max(ncs, 1)
res["rust_time_by_p"] = T["rust"].groupby("p").seconds.median().round(4).to_dict()
if "susier" in T:
    res["susier_time_by_p"] = T["susier"].groupby("p").seconds.median().round(4).to_dict()
s = res.get("susier", {})
res["pass_parity"] = bool(s and s["pip_r"] >= 0.999 and s["cs_identical_share"] >= 0.95)
res["pass_speed_susier"] = bool(s and s["median_time_ratio_vs_rust"] >= 10)
res["pass_speed_pysusie"] = bool(res.get("pysusie_default", {}).get("median_time_ratio_vs_rust", 0) >= 3)
json.dump(res, open(f"{B}/phase1_gate.json", "w"), indent=1)
print(json.dumps(res, indent=1))
