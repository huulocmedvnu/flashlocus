"""GRCh38 variant keys for the phase-1 real loci (fit order = snps.tsv order): chr_pos38_a_b, alleles sorted.
PCOS, FIB, ENDO: GRCh38 position from the source window file by rsid. ANM: hg19 bim positions lifted with
UCSC liftOver. Unmapped and duplicated keys are written as NA (they never match a QTL variant)."""
import os, subprocess, sys
import pandas as pd
S = "/mnt/hdd2/loc-tran/susie_rs"
T = "/mnt/hdd2/loc-tran/repro_coloc/tools"
out = f"{S}/atlas/keys"; os.makedirs(out, exist_ok=True)
rows = []
for trait in ["PCOS", "FIB", "ENDO", "ANM"]:
    win = dict(l.rstrip("\n").split("\t") for l in open(f"{S}/real/win_{trait}.tsv"))
    ids = sorted(win)
    snps = {i: pd.read_csv(f"{S}/real/{trait}/{i}.snps.tsv", sep="\t", dtype={"chr": str}) for i in ids if os.path.exists(f"{S}/real/{trait}/{i}.snps.tsv")}
    if trait == "ANM":
        allv = pd.concat(snps.values()).drop_duplicates(["chr", "pos"])
        bed = f"{out}/anm19.bed"
        with open(bed, "w") as f:
            for c, p in zip(allv.chr, allv.pos): f.write(f"chr{c}\t{p-1}\t{p}\t{c}:{p}\n")
        subprocess.run([f"{T}/liftOver", bed, f"{T}/hg19ToHg38.over.chain.gz", f"{out}/anm38.bed", f"{out}/anm.unmapped"], check=True)
        lo = pd.read_csv(f"{out}/anm38.bed", sep="\t", header=None, names=["c", "s", "e", "k"])
        lo = lo[lo.c.str.replace("chr", "") == lo.k.str.split(":").str[0]]
        m19 = dict(zip(lo.k, lo.e))
    for i, d in snps.items():
        if trait == "ANM":
            p38 = [m19.get(f"{c}:{p}") for c, p in zip(d.chr, d.pos)]
        else:
            w = pd.read_csv(win[i], sep="\t", dtype=str)
            rs = "rsids" if trait == "ENDO" else ("rsid" if "rsid" in w.columns else "rs_id")
            ps = "pos" if trait == "ENDO" else "base_pair_location"
            w[rs] = w[rs].str.split(",").str[0]
            m = dict(zip(w[rs], w[ps].astype(int)))
            p38 = [m.get(r) for r in d.id]
        keys, seen = [], set()
        for c, p, a, b in zip(d.chr, p38, d.a1, d.a2):
            if p is None: keys.append("NA"); continue
            x, y = sorted([a, b]); k = f"{c}_{int(p)}_{x}_{y}"
            keys.append("NA" if k in seen else k); seen.add(k)
        open(f"{out}/{i}.keys", "w").write("\n".join(keys) + "\n")
        rows.append((i, f"{S}/final/real/rs_auto/{i}", f"{out}/{i}.keys", sum(k == "NA" for k in keys), len(keys)))
pd.DataFrame(rows, columns=["id", "fit", "keys", "na", "p"]).to_csv(f"{S}/atlas/loci_summary.tsv", sep="\t", index=False)
with open(f"{S}/atlas/loci.tsv", "w") as f:
    for r in rows: f.write(f"{r[0]}\t{r[1]}\t{r[2]}\n")
print(len(rows), "loci; NA fraction", sum(r[3] for r in rows) / sum(r[4] for r in rows))
