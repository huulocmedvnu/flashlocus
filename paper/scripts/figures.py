"""Data figures 2-5 for the flashLocus GENETICS manuscript. Style as in the MHR organoid paper: plotly_white,
Plotly colours, Arial, black text, no gridlines, lowercase panel titles, no interpretation inside the figure.
Input: paper/figdata (extracted on the benchmark server by extract_figdata.py). Output: paper/figures (PNG, TIFF)."""
import gzip, io, os, re
import numpy as np, pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from PIL import Image

D = os.path.join(os.path.dirname(__file__), "..", "figdata")
O = os.path.join(os.path.dirname(__file__), "..", "figures"); os.makedirs(O, exist_ok=True)
BLACK, GREY = "#111111", "#9a9a9a"
C = {"flashLocus": "#636EFA", "flashLocus, genotype factor": "#19D3F3", "susieR 0.14.2": "#EF553B",
     "susieR 0.16.6": "#AB63FA", "pysusie": "#00CC96", "PLINK 1.9": "#FFA15A", "PLINK 2": "#FF6692"}
W, FS = 900, 13
TOOL = {"rs_sym": "flashLocus", "rs_auto": "flashLocus, genotype factor", "susier014": "susieR 0.14.2",
        "susier016": "susieR 0.16.6", "pysusie": "pysusie"}


def style(fig, h):
    fig.update_layout(template="plotly_white", width=W, height=h, font=dict(family="Arial", size=FS, color=BLACK),
                      margin=dict(l=70, r=20, t=40, b=60), legend=dict(font=dict(size=FS, color=BLACK)))
    fig.update_xaxes(showgrid=False, zeroline=False, showline=True, linecolor=BLACK, ticks="outside", tickcolor=BLACK,
                     tickfont=dict(color=BLACK), title_font=dict(color=BLACK))
    fig.update_yaxes(showgrid=False, zeroline=False, showline=True, linecolor=BLACK, ticks="outside", tickcolor=BLACK,
                     tickfont=dict(color=BLACK), title_font=dict(color=BLACK))
    for a in fig.layout.annotations:
        if a.text and re.match(r"^[a-z]  ", a.text):
            a.update(x=a.x - 0.0, xanchor="center", font=dict(size=FS + 1, color=BLACK))
    return fig


def fmt(v):
    """0.0123 -> 0.012, 6.5e-08 -> 6.5 x 10^-8 written plainly"""
    if v >= 0.01: return f"{v:.2f}"
    m, e = f"{v:.1e}".split("e")
    return f"{m} × 10<sup>{int(e)}</sup>"


LOGT = [1e-3, 1e-2, 1e-1, 1, 10, 100, 1000]
LOGL = ["0.001", "0.01", "0.1", "1", "10", "100", "1000"]


def save(fig, name):
    png = os.path.join(O, f"{name}.png")
    fig.write_image(png, scale=4)
    Image.open(png).convert("RGB").save(os.path.join(O, f"{name}.tif"), compression="tiff_lzw", dpi=(600, 600))
    print("wrote", name)


def fig2():
    """PIP agreement with susieR 0.14.2: rows simulated / real, columns flashLocus, pysusie, susieR 0.16.6."""
    cols = [("rs_sym", "flashLocus"), ("pysusie", "pysusie"), ("susier016", "susieR 0.16.6")]
    titles = [f"{l}  {n}, {s} loci" for (s, row) in [("simulated", 0), ("real", 1)] for l, (_, n) in zip("abc" if row == 0 else "def", cols)]
    fig = make_subplots(rows=2, cols=3, subplot_titles=titles, horizontal_spacing=0.09, vertical_spacing=0.17)
    rng = np.random.default_rng(1)
    for r, s in enumerate(["sim", "real"], 1):
        p = pd.read_csv(f"{D}/pip_{s}.tsv.gz", sep="\t")
        mx = pd.read_csv(f"{D}/pipmaxdiff_{s}.tsv", sep="\t")
        for c, (k, n) in enumerate(cols, 1):
            q = p[p.tool == k]
            far = q[(q.pip_tool - q.pip_susier014).abs() > 1e-3]
            near = q.drop(far.index)
            near = near.iloc[rng.choice(len(near), min(len(near), 6000), replace=False)]
            q = pd.concat([near, far])
            fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode="lines", line=dict(color=GREY, width=1), showlegend=False, hoverinfo="skip"), r, c)
            fig.add_trace(go.Scatter(x=q.pip_susier014, y=q.pip_tool, mode="markers", showlegend=False,
                                     marker=dict(size=3.5, color=C[n], opacity=0.55, line=dict(width=0))), r, c)
            full = p[p.tool == k]
            rr = np.corrcoef(full.pip_susier014, full.pip_tool)[0, 1]
            fig.add_annotation(row=r, col=c, x=0.02, y=0.98, xanchor="left", yanchor="top", showarrow=False, align="left",
                               text=f"max difference {fmt(mx[k].max())}<br>loci > 0.01: {(mx[k] > 0.01).sum()}",
                               font=dict(size=FS - 1, color=BLACK))
            fig.update_xaxes(range=[-0.02, 1.02], tickvals=[0, 0.5, 1], title_text="PIP, susieR 0.14.2" if r == 2 else None, row=r, col=c)
            fig.update_yaxes(range=[-0.02, 1.02], tickvals=[0, 0.5, 1], title_text=f"PIP, {n}", row=r, col=c)
    style(fig, 640)
    save(fig, "Figure2_pip_agreement")


def fig3():
    """Per-locus time against number of variants (one pinned core) and 48-core wall time."""
    fig = make_subplots(rows=1, cols=3, column_widths=[0.35, 0.35, 0.30], horizontal_spacing=0.10,
                        subplot_titles=["a  Simulated loci, one core", "b  Real loci, one core", "c  All loci, 48 cores"])
    order = ["susier014", "susier016", "pysusie", "rs_sym", "rs_auto"]
    for c, s in enumerate(["sim", "real"], 1):
        t = pd.read_csv(f"{D}/timing_{s}.tsv", sep="\t")
        for k in order:
            n = TOOL[k]
            fig.add_trace(go.Scatter(x=t.p, y=t[k], mode="markers", name=n, legendgroup=n, showlegend=(c == 1),
                                     marker=dict(size=4, color=C[n], opacity=0.7, line=dict(width=0))), 1, c)
        fig.update_xaxes(type="log", title_text="Variants in locus", tickvals=[500, 1000, 2000, 5000, 10000],
                         ticktext=["500", "1000", "2000", "5000", "10 000"], row=1, col=c)
        fig.update_yaxes(type="log", title_text="Seconds per locus" if c == 1 else None, tickvals=LOGT, ticktext=LOGL, row=1, col=c)
    tp = {"sim": {"susier014": 71.8, "susier016": 88.8, "pysusie": 16.4, "rs_sym": 4.79, "rs_auto": 1.71},
          "real": {"susier014": 448.1, "susier016": 282.7, "pysusie": 86.8, "rs_sym": 44.1, "rs_auto": 11.9}}
    for k in order:
        n = TOOL[k]
        fig.add_trace(go.Bar(x=["simulated<br>(200)", "real<br>(396)"], y=[tp["sim"][k], tp["real"][k]], name=n, legendgroup=n,
                             showlegend=False, marker_color=C[n]), 1, 3)
    fig.update_yaxes(type="log", title_text="Wall time, seconds", tickvals=LOGT, ticktext=LOGL, row=1, col=3)
    fig.update_layout(barmode="group", legend=dict(orientation="h", y=-0.22, x=0.5, xanchor="center"))
    style(fig, 470)
    save(fig, "Figure3_speed")


def fig4():
    """Colocalization: PP.H4 flashLocus against coloc 5.2.3, and wall time per tissue."""
    p = pd.read_csv(f"{D}/coloc_pairs.tsv.gz", sep="\t")
    names = {"uterus": "uterus", "ovary": "ovary", "vagina": "vagina", "adipose_visceral": "adipose, visceral",
             "blood": "whole blood", "adipose_subcut": "adipose, subcutaneous"}
    fl = pd.read_csv(f"{D}/atlas_flashlocus_timing.tsv", sep="\t").set_index("tissue")
    rt = {}
    for l in open(f"{D}/atlas_r_logs.txt"):
        m = re.match(r"R_(\S+) wall ([\d.]+) s", l)
        if m: rt[m.group(1)] = float(m.group(2))  # last occurrence wins: the rerun of the three large tissues
    fig = make_subplots(rows=1, cols=2, column_widths=[0.42, 0.58], horizontal_spacing=0.12,
                        subplot_titles=["a  PP.H4 per credible-set pair", "b  Wall time per tissue, 48 cores"])
    fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode="lines", line=dict(color=GREY, width=1), showlegend=False), 1, 1)
    fig.add_trace(go.Scatter(x=p["PP.H4_r"], y=p["PP.H4_fl"], mode="markers", showlegend=False,
                             marker=dict(size=4, color=C["flashLocus"], opacity=0.5, line=dict(width=0))), 1, 1)
    d = (p["PP.H4_fl"] - p["PP.H4_r"]).abs().max()
    fig.add_annotation(row=1, col=1, x=0.02, y=0.98, xanchor="left", yanchor="top", showarrow=False, align="left",
                       text=f"{len(p):,} tested pairs<br>max difference {fmt(d)}", font=dict(size=FS - 1, color=BLACK))
    fig.update_xaxes(range=[-0.02, 1.02], tickvals=[0, 0.5, 1], title_text="PP.H4, coloc 5.2.3 (R)", row=1, col=1)
    fig.update_yaxes(range=[-0.02, 1.02], tickvals=[0, 0.5, 1], title_text="PP.H4, flashLocus", row=1, col=1)
    ts = list(names)
    fig.add_trace(go.Bar(x=[names[t] for t in ts], y=[rt[t] for t in ts], name="coloc 5.2.3 (R)", marker_color=C["susieR 0.14.2"]), 1, 2)
    fig.add_trace(go.Bar(x=[names[t] for t in ts], y=[fl.loc[t, "load_s"] + fl.loc[t, "coloc_s"] for t in ts], name="flashLocus",
                         marker_color=C["flashLocus"]), 1, 2)
    fig.update_yaxes(title_text="Seconds", row=1, col=2)
    fig.update_xaxes(tickangle=-30, row=1, col=2)
    fig.update_layout(barmode="group", legend=dict(x=1.0, xanchor="right", y=1.0))
    style(fig, 470)
    save(fig, "Figure4_coloc")


def fig5():
    """LD: wall time per window against PLINK, and largest absolute difference per window."""
    t = pd.read_csv(f"{D}/ld_times.tsv", sep="\t")
    chrx = {"ANM0236", "ANM0237", "ANM0238", "ANM0239", "ANM0240", "ANM0241", "ENDO0036"}
    # chrX after the PLINK 1.9 haploid-male coding (chrx_check.sh, 31/10/2026)
    xfix = {"ANM0236": 2.22e-16, "ANM0237": 3.33e-16, "ANM0238": 2.22e-16, "ANM0239": 3.33e-16, "ANM0240": 2.22e-16,
            "ANM0241": 2.22e-16, "ENDO0036": 2.22e-16}
    fig = make_subplots(rows=1, cols=2, column_widths=[0.6, 0.4], horizontal_spacing=0.13,
                        subplot_titles=["a  Wall time per window, 48 threads", "b  Largest absolute difference in r"])
    for k, n in [("plink19_s", "PLINK 1.9"), ("plink2_s", "PLINK 2"), ("susiers_s", "flashLocus")]:
        fig.add_trace(go.Scatter(x=t.p, y=t[k], mode="markers", name=n, marker=dict(size=4.5, color=C[n], opacity=0.75, line=dict(width=0))), 1, 1)
    fig.update_xaxes(type="log", title_text="Variants in window", tickvals=[500, 1000, 2000, 5000, 10000],
                     ticktext=["500", "1000", "2000", "5000", "10 000"], row=1, col=1)
    fig.update_yaxes(type="log", title_text="Seconds", tickvals=[0.5, 1, 2, 5, 10], ticktext=["0.5", "1", "2", "5", "10"], row=1, col=1)
    a = t[~t.id.isin(chrx)]
    groups = [("PLINK 1.9<br>autosomes", a.maxdiff_p19.astype(float)), ("PLINK 1.9<br>chrX", pd.Series(list(xfix.values()))),
              ("PLINK 2<br>autosomes", a.maxdiff_p2_abs.astype(float))]
    rng = np.random.default_rng(2)
    for i, (lab, v) in enumerate(groups):
        v = v.clip(lower=1e-17)
        fig.add_trace(go.Scatter(x=i + rng.uniform(-0.18, 0.18, len(v)), y=v, mode="markers", showlegend=False,
                                 marker=dict(size=4.5, color=C["flashLocus"], opacity=0.6, line=dict(width=0))), 1, 2)
    fig.update_xaxes(tickvals=[0, 1, 2], ticktext=[g[0] for g in groups], range=[-0.6, 2.6], row=1, col=2)
    fig.update_yaxes(type="log", range=[-17, -14], tickvals=[1e-17, 1e-16, 1e-15, 1e-14],
                     ticktext=["10<sup>-17</sup>", "10<sup>-16</sup>", "10<sup>-15</sup>", "10<sup>-14</sup>"], title_text="Absolute difference", row=1, col=2)
    fig.update_layout(legend=dict(x=0.02, y=1.0))
    style(fig, 450)
    save(fig, "Figure5_ld")


if __name__ == "__main__":
    fig2(); fig3(); fig4(); fig5()
