"""Throughput: pysusie parallel over loci (multiprocessing, one Numba thread per worker), wall time incl. I/O."""
import sys, time, os
os.environ.setdefault("NUMBA_NUM_THREADS", "1")
import numpy as np
from multiprocessing import Pool


def init():
    global SuSiE
    from pysusie import SuSiE
    z = np.zeros(50); z[3] = 6.0; R = np.eye(50)
    SuSiE(n_effects=10).fit_from_summary_stats(z=z, R=R, n=10000.0)   # JIT warm-up per worker


def one(row):
    lid, zf, lf, n = row
    z = np.loadtxt(zf); p = len(z); R = np.fromfile(lf, dtype="<f8").reshape(p, p)
    try:
        m = SuSiE(n_effects=10).fit_from_summary_stats(z=z, R=R, n=float(n))
        return len(m.result_.get_credible_sets(R=R))
    except Exception:
        return -1


if __name__ == "__main__":
    man, cores = sys.argv[1], int(sys.argv[2])
    rows = [l.rstrip("\n").split("\t") for l in open(man) if l.strip()]
    with Pool(cores, initializer=init) as pool:
        pool.map(one, rows[:cores])                                      # make sure every worker is warm
        t0 = time.perf_counter()
        pool.map(one, rows, chunksize=1)
        print("pysusie loci", len(rows), "cores", cores, "wall_s", time.perf_counter() - t0)
