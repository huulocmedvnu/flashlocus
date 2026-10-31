# Data and scripts for the manuscript figures

- `figdata/`: per-variant PIPs, per-locus timings, coloc pairs, LD timings and throughput logs, extracted
  from the benchmark outputs by `scripts/extract_figdata.py` (set `BENCH_ROOT` to the benchmark directory).
- `scripts/figures.py`: Figures 2 to 5 from `figdata/` (Python with plotly and kaleido).
- `scripts/figure1_schematic.py`: Figure 1 as SVG (needs `rsvg-convert`).
- `Table_S1_colocalization_pairs.csv`: the 140 credible-set pairs with PP.H4 of at least 0.8.

    python paper/scripts/figures.py && python paper/scripts/figure1_schematic.py   # writes paper/figures/
