"""Phase 1 benchmark loci (PLAN.md): 200 simulated loci from 1000 Genomes EUR (503 samples).
Per locus: contiguous block of p variants with MAF >= 1%, R = in-sample correlation, k in {1,2,3} causal
variants, total variance explained h2 ~ U(0.0005, 0.002) split evenly, GWAS n = 100,000.
z = sqrt(n) * R[:, C] @ (s * sqrt(h2/k)) + X' g / sqrt(N), g ~ N(0, I_N), so Cov(noise) = R exactly."""
import os, sys, json
import numpy as np, pandas as pd

REF = "/mnt/hdd2/loc-tran/infertility_work/ref/g1000_eur"
OUT = "/mnt/hdd2/loc-tran/susie_rs/bench/sim"
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(20261031)
N = 503
NB = (N + 3) // 4
bim = pd.read_csv(REF + ".bim", sep="\t", header=None, usecols=[0, 1, 3], names=["chr", "snp", "pos"], dtype={"chr": str})
bed = np.memmap(REF + ".bed", dtype=np.uint8, mode="r", offset=3).reshape(len(bim), NB)
LUT = np.array([2.0, np.nan, 1.0, 0.0])  # PLINK: 00 hom A1, 01 missing, 10 het, 11 hom A2 -> count of A1


def geno(i0, i1):
    b = np.asarray(bed[i0:i1])
    codes = np.stack([(b >> (2 * k)) & 3 for k in range(4)], axis=2).reshape(i1 - i0, NB * 4)[:, :N]
    g = LUT[codes]
    m = np.nanmean(g, axis=1, keepdims=True)
    g = np.where(np.isnan(g), m, g)
    return g


sizes = [500, 1000, 2000, 3000, 5000]
auto = bim.index[bim.chr.isin([str(c) for c in range(1, 23)])].values
rows = []
L = 0
while L < 200:
    p = sizes[L % len(sizes)]
    s = int(rng.choice(auto[: -p * 4]))
    g = geno(s, s + p * 4)
    if bim.chr.iloc[s] != bim.chr.iloc[s + p * 4 - 1]:
        continue
    af = g.mean(1) / 2
    keep = np.flatnonzero((np.minimum(af, 1 - af) >= 0.01) & (g.std(1) > 0))
    if len(keep) < p:
        continue
    keep = keep[:p]
    X = g[keep].T
    X = (X - X.mean(0)) / X.std(0)
    R = (X.T @ X) / N
    np.fill_diagonal(R, 1.0)
    k = 1 + L % 3
    h2 = rng.uniform(0.0005, 0.002)
    C = rng.choice(p, size=k, replace=False)
    beta = rng.choice([-1, 1], size=k) * np.sqrt(h2 / k)
    n = 100000
    z = np.sqrt(n) * R[:, C] @ beta + X.T @ rng.standard_normal(N) / np.sqrt(N)
    lid = f"sim{L:03d}"
    np.savetxt(f"{OUT}/{lid}.z", z, fmt="%.8f")
    R.astype("<f8").tofile(f"{OUT}/{lid}.ld.bin")
    rows.append({"id": lid, "p": p, "k": k, "h2": h2, "causal": ",".join(str(c + 1) for c in sorted(C)), "chr": bim.chr.iloc[s],
                 "start": int(bim.pos.iloc[s + keep[0]]), "n": n, "max_abs_z": float(np.abs(z).max())})
    L += 1
    print(lid, p, k, round(float(np.abs(z).max()), 1), flush=True)
m = pd.DataFrame(rows)
m.to_csv(f"{OUT}/truth.tsv", sep="\t", index=False)
with open(f"{OUT}/manifest.tsv", "w") as f:
    for r in rows:
        f.write(f"{r['id']}\t{OUT}/{r['id']}.z\t{OUT}/{r['id']}.ld.bin\t{r['n']}\n")
print("DONE", len(rows))
