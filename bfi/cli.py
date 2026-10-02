"""
BFI command-line interface.

Usage:
    bfi run --bfile PREFIX --ne 911 --recomb-rate 4e-8 --keep indica.txt --out bfi.tsv
    bfi run --bfile PREFIX --ne 1263 --recomb-rate 2.5e-8 --reference improved.txt --out cross.tsv
"""

import argparse
import sys

from bfi import __version__
from bfi.core import compute_bfi, EQ_MODES


def main():
    parser = argparse.ArgumentParser(prog="bfi", description="Breeding Freedom Index")
    parser.add_argument("--version", action="version", version=f"bfi {__version__}")
    sub = parser.add_subparsers(dest="command")

    r = sub.add_parser("run", help="Compute BFI")
    r.add_argument("--bfile", required=True, help="PLINK binary prefix (.bed/.bim/.fam), unique variant IDs")
    r.add_argument("--ne", type=float, required=True, help="Effective population size of the reference group")
    r.add_argument("--recomb-rate", type=float, default=4e-8, help="Per-bp recombination rate (default 4e-8)")
    r.add_argument("--keep", default=None, help="Samples to score (IID or FID IID per line); default all")
    r.add_argument("--reference", default=None,
                   help="Reference samples defining blocks/major alleles/weights (cross-group BFI); "
                        "default: the scored samples")
    r.add_argument("--chromosomes", default=None, help="Comma-separated chromosomes (default: all)")
    r.add_argument("--blocks-max-kb", type=int, default=200, help="PLINK --blocks-max-kb (default 200)")
    r.add_argument("--blocks-min-maf", type=float, default=0.05, help="PLINK --blocks-min-maf (default 0.05)")
    r.add_argument("--maf", type=float, default=None, help="Optional --maf filter used with --blocks")
    r.add_argument("--eq-mode", choices=EQ_MODES, default="span",
                   help="Expected r2: 'span' = equilibrium r2 at the block span BP2-BP1 (default, "
                        "published definition); 'pair' = averaged over the SNP pairs of r2_obs (sensitivity)")
    r.add_argument("--r2-mode", choices=("pairwise", "complete"), default="pairwise",
                   help="Missing data in r2_obs: 'pairwise' = per SNP pair, accessions genotyped at both "
                        "(default, v1.0.2); 'complete' = only SNPs without missing calls (v1.0.0-1.0.1)")
    r.add_argument("--sampling-term", action="store_true", help="Add 1/n to the expected r2")
    r.add_argument("--workers", type=int, default=6, help="Parallel workers (default 6)")
    r.add_argument("--plink", default="plink", help="PLINK 1.9 executable")
    r.add_argument("--blocks-out", default=None, help="Also write the block table (TSV)")
    r.add_argument("--out", default="bfi_results.tsv", help="Output TSV (IID, BFI)")
    a = parser.parse_args()

    if a.command is None:
        parser.print_help()
        sys.exit(1)

    chroms = a.chromosomes.split(",") if a.chromosomes else None
    print(f"bfi {__version__} | {a.bfile} | Ne={a.ne} | c={a.recomb_rate}/bp | eq-mode={a.eq_mode}"
          f"{' +1/n' if a.sampling_term else ''} | r2-mode={a.r2_mode} | blocks-max-kb={a.blocks_max_kb}")
    res = compute_bfi(bfile=a.bfile, sample_file=a.keep, reference_file=a.reference, ne=a.ne,
                      recomb_rate=a.recomb_rate, chromosomes=chroms, blocks_max_kb=a.blocks_max_kb,
                      blocks_min_maf=a.blocks_min_maf, maf=a.maf, eq_mode=a.eq_mode, r2_mode=a.r2_mode,
                      sampling_term=a.sampling_term, n_workers=a.workers, plink_path=a.plink,
                      return_blocks=bool(a.blocks_out))
    if a.blocks_out:
        res, blocks = res
        blocks.to_csv(a.blocks_out, sep="\t", index=False, float_format="%.6f")
    res.to_csv(a.out, sep="\t", index=False, float_format="%.6f")
    print(f"  {len(res)} samples -> {a.out} | BFI {res.BFI.min():.4f}-{res.BFI.max():.4f}, mean {res.BFI.mean():.4f}")


if __name__ == "__main__":
    main()
