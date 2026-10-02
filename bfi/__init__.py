"""
BFI - Breeding Freedom Index

A quantitative metric for genomic design space in crop breeding.
Measures an individual's potential to disrupt excess LD blocks
when used as a crossing parent.

Reference:
    Kim, M., Lee, K., Kim, K.D., Hwang, J.-H., & Lee, T.-H. (2026).
    Breeding Freedom Index: an excess-linkage-disequilibrium-weighted framework
    for genomic design space in crop breeding.
"""

__version__ = "1.0.2"

from bfi.core import (compute_bfi, compute_bfi_chromosome, block_r2_obs, block_r2_eq,
                      major_homozygote, rare_allele_burden)
from bfi.ne import estimate_ne, equilibrium_r2
from bfi.ld import define_blocks, compute_ld_decay
