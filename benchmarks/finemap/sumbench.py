import os, numpy as np, pandas as pd
def cs_r(f):
    if not os.path.exists(f) or os.path.getsize(f)==0: return set()
    return {tuple(sorted(map(int,x[1].split(",")))) for x in (l.rstrip("\n").split("\t") for l in open(f)) if len(x)>1}
def cs_rs(f):
    L=[l.rstrip("\n").split("\t") for l in open(f)][1:]
    return {tuple(sorted(map(int,x[6].split(",")))) for x in L if len(x)>6}
for s,man in [("sim","../bench/sim2/manifest.tsv"),("real","../real/all.tsv")]:
    ids=[l.split("\t")[0] for l in open(man) if l.strip()]
    T={t:pd.read_csv(f"{s}/{t}/timing.tsv",sep="\t").set_index("id").seconds for t in ["susier014","susier016","pysusie","pysusie_exact","rs_sym","rs_auto"]}
    print("#####",s,"loci",len(ids), "| susieR 0.14.2 total_s %.1f median_s %.3f" % (T["susier014"].sum(), T["susier014"].median()))
    for t in ["susier016","pysusie","pysusie_exact","rs_sym","rs_auto"]:
        a=[];b=[];same=0;big=0
        for i in ids:
            x=np.loadtxt(f"{s}/{t}/{i}.pip"); y=np.loadtxt(f"{s}/susier014/{i}.pip")
            a.append(x); b.append(y); big+= np.abs(x-y).max()>0.01
            cs=cs_rs(f"{s}/{t}/{i}.cs.tsv") if t.startswith("rs") else cs_r(f"{s}/{t}/{i}.cs")
            same+= cs==cs_r(f"{s}/susier014/{i}.cs")
        a=np.concatenate(a); b=np.concatenate(b)
        r=(T[t]/T["susier014"]).median()
        print(f"{t:14s} pip_r={np.corrcoef(a,b)[0,1]:.4f} loci_pipdiff>0.01={big:3d} cs_identical={same/len(ids):.3f} total_s={T[t].sum():8.1f} median_speedup_vs_susieR014={1/r:6.1f}")
