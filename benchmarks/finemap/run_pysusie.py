"""pysusie on the same loci. mode 'default' = package defaults (prior variance optimised every 5th iteration);
mode 'exact' = optimised every iteration (prior_variance_optim_period=1), the susieR algorithm."""
import sys, os, time
import numpy as np, pandas as pd
from pysusie import SuSiE

man, out, mode = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(out, exist_ok=True)
rows = [l.rstrip("\n").split("\t") for l in open(man) if l.strip()]
import json
kw = dict(n_effects=10)
if mode == "exact":
    kw["prior_variance_optim_period"] = 1
elif mode.startswith("{"):          # explicit configuration, e.g. '{"xtx_update_mode": "batched", "xtx_batch_size": 4}'
    kw.update(json.loads(mode))
# JIT warm-up on the first locus (not timed)
z0 = np.loadtxt(rows[0][1]); R0 = np.fromfile(rows[0][2], dtype="<f8").reshape(len(z0), len(z0))
SuSiE(**kw).fit_from_summary_stats(z=z0, R=R0, n=float(rows[0][3]))
tim = []
for lid, zf, lf, n in rows:
    z = np.loadtxt(zf); p = len(z); R = np.fromfile(lf, dtype="<f8").reshape(p, p)
    t = time.perf_counter()
    try:
        m = SuSiE(**kw).fit_from_summary_stats(z=z, R=R, n=float(n))
        res = m.result_
        sets = res.get_credible_sets(R=R)   # purity filtering timed, as in susieR
        pip = np.asarray(res.pip)
        s = time.perf_counter() - t
        np.savetxt(f"{out}/{lid}.pip", pip, fmt="%.10e")
        with open(f"{out}/{lid}.cs", "w") as f:
            for c in sets:
                f.write(f"{int(c.effect_index) + 1}\t{','.join(str(int(x) + 1) for x in sorted(c.variables))}\n")
        tim.append((lid, p, s, res.n_iter, len(sets), ""))
    except Exception as e:
        tim.append((lid, p, np.nan, -1, -1, repr(e)[:200]))
pd.DataFrame(tim, columns=["id", "p", "seconds", "niter", "n_cs", "err"]).to_csv(f"{out}/timing.tsv", sep="\t", index=False)
