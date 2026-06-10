"""
BFI core computation engine.

BFI = Σ_b [ RAB_b × w_b ] / Σ_b [ w_b ]

where:
    RAB_b = rare allele burden in block b (non-conformity to majority allele)
    w_b   = max(r²_obs,b - r²_eq,b, 0)  (excess LD weight)
"""

import os
import subprocess
import shutil
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

from bfi.ne import equilibrium_r2


def compute_bfi(
    bfile,
    sample_file=None,
    ne=None,
    recomb_rate=4e-8,
    chromosomes=None,
    blocks_max_kb=200,
    blocks_min_maf=0.05,
    min_block_snps=2,
    n_workers=6,
    plink_path="plink",
    tmpdir=None,
):
    """
    Compute genome-wide BFI for all individuals in a PLINK binary fileset.

    Parameters
    ----------
    bfile : str
        Path prefix for PLINK binary files (.bed/.bim/.fam).
    sample_file : str, optional
        File with sample IDs to include (one per line, or FID\\tIID format).
    ne : float
        Effective population size. Required.
    recomb_rate : float
        Per-base-pair recombination rate. Default: 4e-8 (rice).
    chromosomes : list of int, optional
        Chromosomes to process. Default: auto-detect from .bim file.
    blocks_max_kb : int
        Maximum LD block size in kb for PLINK --blocks.
    blocks_min_maf : float
        Minimum MAF for block SNPs.
    min_block_snps : int
        Minimum SNPs per block.
    n_workers : int
        Number of parallel workers (one per chromosome).
    plink_path : str
        Path to PLINK 1.9 executable.
    tmpdir : str, optional
        Temporary directory. Default: creates one next to bfile.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ['IID', 'BFI'] sorted by BFI descending.
    """
    if ne is None:
        raise ValueError("ne (effective population size) is required")

    if chromosomes is None:
        bim = pd.read_csv(f"{bfile}.bim", sep=r'\s+', header=None,
                          names=["CHR", "SNP", "CM", "BP", "A1", "A2"])
        chromosomes = sorted(bim['CHR'].unique())

    if tmpdir is None:
        tmpdir = os.path.join(os.path.dirname(bfile), "_bfi_tmp")

    args_list = [
        (bfile, sample_file, ne, recomb_rate, chr_num,
         blocks_max_kb, blocks_min_maf, min_block_snps, plink_path, tmpdir)
        for chr_num in chromosomes
    ]

    chr_results = []
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(compute_bfi_chromosome, *args): args[4]
                   for args in args_list}
        for future in as_completed(futures):
            chr_num = futures[future]
            result = future.result()
            if result is not None:
                chr_results.append(result)

    if not chr_results:
        raise RuntimeError("No blocks found on any chromosome")

    # Merge across chromosomes
    merged = chr_results[0][['IID']].copy()
    merged['w_bfi'] = sum(r['w_bfi'].values for r in chr_results)
    merged['w_total'] = sum(r['w_total'].values for r in chr_results)
    merged['BFI'] = np.where(merged['w_total'] > 0,
                              merged['w_bfi'] / merged['w_total'], 0.0)

    return merged[['IID', 'BFI']].sort_values('BFI', ascending=False).reset_index(drop=True)


def compute_bfi_chromosome(
    bfile, sample_file, ne, recomb_rate, chr_num,
    blocks_max_kb, blocks_min_maf, min_block_snps, plink_path, tmpdir_base,
):
    """
    Compute per-chromosome BFI components (w_bfi, w_total) for all individuals.

    Returns DataFrame with ['IID', 'w_bfi', 'w_total'] or None if no blocks found.
    """
    tmpdir = os.path.join(tmpdir_base, f"chr{chr_num}")
    os.makedirs(tmpdir, exist_ok=True)

    try:
        return _compute_chr_inner(
            bfile, sample_file, ne, recomb_rate, chr_num,
            blocks_max_kb, blocks_min_maf, min_block_snps, plink_path, tmpdir
        )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def _compute_chr_inner(
    bfile, sample_file, ne, recomb_rate, chr_num,
    blocks_max_kb, blocks_min_maf, min_block_snps, plink_path, tmpdir,
):
    # Keep file
    keep_args = []
    if sample_file:
        keep_file = os.path.join(tmpdir, "keep.txt")
        with open(sample_file) as fin:
            samples = [line.strip() for line in fin if line.strip()]
        with open(keep_file, "w") as f:
            for s in samples:
                if "\t" in s:
                    f.write(s + "\n")
                else:
                    f.write(f"0\t{s}\n")
        keep_args = ["--keep", keep_file]

    # Step 1: LD block definition
    block_prefix = os.path.join(tmpdir, "blocks")
    cmd = [plink_path, "--bfile", bfile] + keep_args + [
        "--chr", str(chr_num), "--blocks", "no-pheno-req",
        "--blocks-max-kb", str(blocks_max_kb),
        "--blocks-min-maf", str(blocks_min_maf),
        "--out", block_prefix,
    ]
    subprocess.run(cmd, capture_output=True, text=True)

    block_file = f"{block_prefix}.blocks.det"
    if not os.path.exists(block_file) or os.path.getsize(block_file) == 0:
        return None

    blocks_df = pd.read_csv(block_file, sep=r'\s+')
    if len(blocks_df) == 0:
        return None

    # Collect SNPs
    all_snps = set()
    block_info = []
    for _, row in blocks_df.iterrows():
        snps = row["SNPS"].split("|")
        block_info.append((snps, row["BP1"], row["BP2"]))
        all_snps.update(snps)

    # Step 2: Extract genotypes
    snp_file = os.path.join(tmpdir, "snps.txt")
    with open(snp_file, "w") as f:
        f.write("\n".join(sorted(all_snps)) + "\n")

    geno_prefix = os.path.join(tmpdir, "geno")
    cmd = [plink_path, "--bfile", bfile] + keep_args + [
        "--extract", snp_file, "--chr", str(chr_num),
        "--recode", "A", "--out", geno_prefix,
    ]
    subprocess.run(cmd, capture_output=True, text=True)

    raw_file = f"{geno_prefix}.raw"
    if not os.path.exists(raw_file):
        return None

    geno = pd.read_csv(raw_file, sep=r'\s+')
    n_samples = len(geno)
    geno_cols = [c for c in geno.columns
                 if c not in ("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")]
    col_to_snp = {"_".join(c.split("_")[:-1]): c for c in geno_cols}

    # Step 3-6: Per-block BFI computation
    w_bfi = np.zeros(n_samples)
    w_total = 0.0

    for snps, bp1, bp2 in block_info:
        cols = [col_to_snp[s] for s in snps if s in col_to_snp]
        if len(cols) < min_block_snps:
            continue

        mat = geno[cols].values
        is_valid = ~np.isnan(mat)
        valid_all = np.all(is_valid, axis=0)
        mat_clean = mat[:, valid_all]
        if mat_clean.shape[1] < 2:
            continue

        # Step 3: Observed r² per block
        stds = np.nanstd(mat_clean, axis=0)
        good = stds > 0
        if good.sum() < 2:
            continue
        mat_std = (mat_clean[:, good] - np.nanmean(mat_clean[:, good], axis=0)) / stds[good]
        corr = np.corrcoef(mat_std.T)
        triu = np.triu_indices(corr.shape[0], k=1)
        obs_r2 = np.mean(corr[triu] ** 2)

        # Step 4: Excess LD weight
        eq_r2 = equilibrium_r2(bp2 - bp1, ne, recomb_rate)
        excess = max(obs_r2 - eq_r2, 0.0)
        if excess <= 0:
            continue

        # Step 5: Rare allele burden (non-conformity)
        snp_means = np.nanmean(mat, axis=0)
        major = np.where(snp_means <= 1.0, 0, 2)
        matches = (mat == major[np.newaxis, :]) & is_valid
        valid_count = np.sum(is_valid, axis=1)
        non_conformity = np.where(
            valid_count > 0,
            1.0 - np.sum(matches, axis=1) / valid_count,
            0.0,
        )

        # Step 6: Accumulate
        w_bfi += non_conformity * excess
        w_total += excess

    return pd.DataFrame({
        "IID": geno["IID"],
        "w_bfi": w_bfi,
        "w_total": w_total,
    })
