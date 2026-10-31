"""flashLocus: summary-statistic fine-mapping (exact port of susieR 0.14.2 susie_rss),
colocalisation (exact port of coloc 5.2.3 coloc.bf_bf) and exact LD from PLINK files, in Rust."""
import sys

from ._core import __version__, coloc_bf_bf, ld, susie_rss
from ._core import cli as _cli

__all__ = ["susie_rss", "coloc_bf_bf", "ld", "main", "__version__"]


def main():
    """Console entry point: `flashlocus <command> ...`."""
    _cli(["flashlocus"] + sys.argv[1:])
