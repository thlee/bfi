"""
BFI command-line interface.

Usage:
    bfi run --bfile PREFIX --ne 911 --recomb-rate 4e-8 --out results.tsv
    bfi run --bfile PREFIX --ne 911 --chromosomes 1,2,3 --workers 12
"""

import argparse
import sys

from bfi.core import compute_bfi


def main():
    parser = argparse.ArgumentParser(
        prog="bfi",
        description="Breeding Freedom Index — genomic design space metric",
    )
    subparsers = parser.add_subparsers(dest="command")

    # run command
    run_parser = subparsers.add_parser("run", help="Compute BFI")
    run_parser.add_argument("--bfile", required=True,
                            help="PLINK binary file prefix (.bed/.bim/.fam)")
    run_parser.add_argument("--ne", type=float, required=True,
                            help="Effective population size")
    run_parser.add_argument("--recomb-rate", type=float, default=4e-8,
                            help="Per-bp recombination rate (default: 4e-8 for rice)")
    run_parser.add_argument("--keep", default=None,
                            help="File with sample IDs to include")
    run_parser.add_argument("--chromosomes", default=None,
                            help="Comma-separated chromosome list (default: all)")
    run_parser.add_argument("--blocks-max-kb", type=int, default=200,
                            help="Max LD block size in kb (default: 200)")
    run_parser.add_argument("--blocks-min-maf", type=float, default=0.05,
                            help="Min MAF for block SNPs (default: 0.05)")
    run_parser.add_argument("--workers", type=int, default=6,
                            help="Parallel workers (default: 6)")
    run_parser.add_argument("--plink", default="plink",
                            help="Path to PLINK 1.9 executable")
    run_parser.add_argument("--out", default="bfi_results.tsv",
                            help="Output file (default: bfi_results.tsv)")

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    if args.command == "run":
        chromosomes = None
        if args.chromosomes:
            chromosomes = [int(c) for c in args.chromosomes.split(",")]

        print(f"BFI v1.0.0")
        print(f"  Input: {args.bfile}")
        print(f"  N_e: {args.ne}")
        print(f"  Recombination rate: {args.recomb_rate}")
        print(f"  Block max kb: {args.blocks_max_kb}")
        print(f"  Workers: {args.workers}")
        print()

        result = compute_bfi(
            bfile=args.bfile,
            sample_file=args.keep,
            ne=args.ne,
            recomb_rate=args.recomb_rate,
            chromosomes=chromosomes,
            blocks_max_kb=args.blocks_max_kb,
            blocks_min_maf=args.blocks_min_maf,
            n_workers=args.workers,
            plink_path=args.plink,
        )

        result.to_csv(args.out, sep="\t", index=False, float_format="%.6f")
        print(f"  Output: {args.out}")
        print(f"  Samples: {len(result)}")
        print(f"  BFI range: {result['BFI'].min():.4f} - {result['BFI'].max():.4f}")
        print(f"  BFI mean: {result['BFI'].mean():.4f}")


if __name__ == "__main__":
    main()
