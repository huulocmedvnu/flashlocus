"""Figure 1: flashLocus schematic as an editable SVG, exported to PNG and 600 dpi TIFF.
Canvas units 10 per mm, 180 x 86 mm, 8 pt = 28 units. Style follows the HR graphical abstract: thin column dividers,
bold column titles, simple icons, coloured labels, black text, values only."""
import os, subprocess
from PIL import Image

O = os.path.join(os.path.dirname(__file__), "..", "figures"); os.makedirs(O, exist_ok=True)
W, H = 1800, 860
BLUE, RED, GREEN, PURPLE, ORANGE = "#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A"
INK, GREY = "#111111", "#9a9a9a"
FONT = "Helvetica, Arial, sans-serif"
svg = []
add = svg.append


def text(x, y, s, size=28, weight="normal", anchor="start", fill=INK):
    add(f"<text x='{x}' y='{y}' font-family='{FONT}' font-size='{size}' font-weight='{weight}' text-anchor='{anchor}' fill='{fill}'>{s}</text>")


def sup(m, e):
    return f"{m} × 10<tspan dy='-11' font-size='21'>−{e}</tspan><tspan dy='11'> </tspan>"


def arrow(x1, y1, x2, y2, colour=GREY):
    add(f"<line x1='{x1}' y1='{y1}' x2='{x2}' y2='{y2}' stroke='{colour}' stroke-width='4' marker-end='url(#ah)'/>")


add(f"<svg xmlns='http://www.w3.org/2000/svg' width='{W}' height='{H}' viewBox='0 0 {W} {H}'>")
add(f"<defs><marker id='ah' markerWidth='10' markerHeight='10' refX='8' refY='5' orient='auto'>"
    f"<path d='M0,0 L10,5 L0,10 Z' fill='{GREY}'/></marker></defs>")
add(f"<rect width='{W}' height='{H}' fill='white'/>")

cols = [(20, 620, "Inputs"), (630, 1170, "flashLocus"), (1180, 1790, "Outputs")]
for x0, x1, t in cols:
    text((x0 + x1) / 2, 60, t, size=34, weight="bold", anchor="middle")
for x in (625, 1175):
    add(f"<line x1='{x}' y1='30' x2='{x}' y2='{H - 95}' stroke='{GREY}' stroke-width='1.5'/>")

rows = [205, 455, 705]

# ---- inputs
# reference genotypes: grid of 0/1/2 calls
gx, gy = 40, rows[0] - 60
shade = [[0, 1, 2, 1, 0, 0], [1, 2, 2, 1, 0, 1], [0, 0, 1, 2, 2, 1], [2, 1, 0, 0, 1, 2]]
for i, r in enumerate(shade):
    for j, v in enumerate(r):
        add(f"<rect x='{gx + j * 24}' y='{gy + i * 24}' width='21' height='21' fill='{ORANGE}' fill-opacity='{0.15 + 0.35 * v}' stroke='{ORANGE}' stroke-width='1.5'/>")
text(210, rows[0] - 22, "Reference genotypes", weight="bold")
text(210, rows[0] + 14, "PLINK .bed / .bim / .fam", size=26)
text(210, rows[0] + 48, "1000 Genomes, European", size=26)
# GWAS: small Manhattan peak
import math
for k in range(26):
    x = 40 + k * 5.6
    peak = 95 * math.exp(-((k - 15) / 3.2) ** 2)
    y = rows[1] + 45 - peak - (k * 7 % 13)
    add(f"<circle cx='{x:.1f}' cy='{y:.1f}' r='4.5' fill='{RED}' fill-opacity='0.75'/>")
add(f"<line x1='34' y1='{rows[1] + 52}' x2='188' y2='{rows[1] + 52}' stroke='{INK}' stroke-width='2'/>")
text(210, rows[1] - 22, "GWAS summary statistics", weight="bold")
text(210, rows[1] + 14, "z scores, sample size", size=26)
text(210, rows[1] + 48, "396 loci, 4 traits", size=26)
# eQTL Catalogue: stacked tissue cards
for k in range(3):
    add(f"<rect x='{44 + k * 14}' y='{rows[2] - 52 + k * 14}' width='110' height='76' rx='8' fill='white' stroke='{GREEN}' stroke-width='3'/>")
for k in range(4):
    add(f"<line x1='{80}' y1='{rows[2] + 2 + k * 12}' x2='{150}' y2='{rows[2] + 2 + k * 12}' stroke='{GREEN}' stroke-width='3' stroke-opacity='0.6'/>")
text(210, rows[2] - 22, "eQTL Catalogue SuSiE fits", weight="bold")
text(210, rows[2] + 14, "log Bayes factors per gene", size=26)
text(210, rows[2] + 48, "6 GTEx tissues", size=26)

# ---- flashLocus modules
mods = [(ORANGE, "Linkage disequilibrium", "exact integer correlation"),
        (RED, "Fine-mapping", "SuSiE-RSS, as susieR 0.14.2"),
        (GREEN, "Colocalization", "credible-set pairs, as coloc 5.2.3")]
for (c, t1, t2), y in zip(mods, rows):
    add(f"<rect x='675' y='{y - 72}' width='450' height='130' rx='14' fill='{c}' fill-opacity='0.10' stroke='{c}' stroke-width='3.5'/>")
    text(900, y - 18, t1, size=30, weight="bold", anchor="middle")
    text(900, y + 22, t2, anchor="middle")
for y in rows:
    arrow(612, y, 666, y)
arrow(900, rows[0] + 60, 900, rows[1] - 80)
arrow(900, rows[1] + 60, 900, rows[2] - 80)
text(900, H - 34, "Rust, AVX-512 / AVX2 / NEON, Python and command line", size=26, anchor="middle")

# ---- outputs
# LD heat map
lx, ly = 1200, rows[0] - 66
vals = [[1, .8, .5, .2, .1], [.8, 1, .6, .3, .1], [.5, .6, 1, .7, .4], [.2, .3, .7, 1, .8], [.1, .1, .4, .8, 1]]
for i, r in enumerate(vals):
    for j, v in enumerate(r):
        add(f"<rect x='{lx + j * 25}' y='{ly + i * 25}' width='24' height='24' fill='{ORANGE}' fill-opacity='{v:.2f}'/>")
text(1345, rows[0] - 22, "Correlation matrix", weight="bold")
text(1345, rows[0] + 14, "difference to PLINK ≤ " + sup("4.4", 16), size=26)
text(1345, rows[0] + 48, "3.7 × faster than PLINK 2", size=26)
# PIP bars
pips = [.02, .03, .05, .9, .08, .03, .02, .4, .35, .02]
for k, v in enumerate(pips):
    hgt = 100 * v
    add(f"<rect x='{1200 + k * 12}' y='{rows[1] + 45 - hgt:.1f}' width='9' height='{hgt:.1f}' fill='{RED}' fill-opacity='0.8'/>")
add(f"<line x1='1196' y1='{rows[1] + 46}' x2='1324' y2='{rows[1] + 46}' stroke='{INK}' stroke-width='2'/>")
text(1345, rows[1] - 22, "PIP and credible sets", weight="bold")
text(1345, rows[1] + 14, "difference to susieR ≤ " + sup("3.0", 5), size=26)
text(1345, rows[1] + 48, "23 × faster than susieR", size=26)
# PP.H4 tiles
tiles = [[.05, .9, .1], [.2, .02, .85], [.95, .1, .3]]
for i, r in enumerate(tiles):
    for j, v in enumerate(r):
        add(f"<rect x='{1205 + j * 36}' y='{rows[2] - 58 + i * 36}' width='32' height='32' rx='4' fill='{GREEN}' fill-opacity='{0.1 + 0.85 * v:.2f}'/>")
text(1345, rows[2] - 22, "PP.H4 per gene and tissue", weight="bold")
text(1345, rows[2] + 14, "difference to coloc ≤ " + sup("5.0", 12), size=26)
text(1345, rows[2] + 48, "7.8 × faster than coloc in R", size=26)
for y in rows:
    arrow(1134, y, 1192, y)
add("</svg>")

svgp = os.path.join(O, "Figure1_schematic.svg")
open(svgp, "w").write("\n".join(svg))
png = os.path.join(O, "Figure1_schematic.png")
subprocess.run(["rsvg-convert", "-z", "4", "-o", png, svgp], check=True)
Image.open(png).convert("RGB").save(os.path.join(O, "Figure1_schematic.tif"), compression="tiff_lzw", dpi=(600, 600))
print("wrote Figure1_schematic")
