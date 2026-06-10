"""
BFI - Breeding Freedom Index

A quantitative metric for genomic design space in crop breeding.
Measures an individual's potential to disrupt excess LD blocks
when used as a crossing parent.

Reference:
    Lee, T.-H. (2026). Breeding Freedom Index: A Quantitative Metric
    for Genomic Design Space in Crop Breeding.
"""

__version__ = "1.0.0"

from bfi.core import compute_bfi, compute_bfi_chromosome
from bfi.ne import estimate_ne, equilibrium_r2
from bfi.ld import define_blocks, compute_ld_decay
