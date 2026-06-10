"""
LD block definition and LD decay computation.

Wraps PLINK 1.9 for LD calculations.
"""

import os
import subprocess
import numpy as np
import pandas as pd


def define_blocks(
    bfile,
    chr_num,
    plink_path="plink",
    sample_file=None,
    blocks_max_kb=200,
    blocks_min_maf=0.05,
    outdir=None,
):
    """
    Define LD blocks using PLINK --blocks (Gabriel et al., 2002).

    Parameters
    ----------
    bfile : str
        PLINK binary file prefix.
    chr_num : int
        Chromosome number.
    plink_path : str
        Path to PLINK 1.9.
    sample_file : str, optional
        Keep file for subsetting samples.
    blocks_max_kb : int
        Maximum block size in kb.
    blocks_min_maf : float
        Minimum MAF for block SNPs.
    outdir : str, optional
        Output directory.

    Returns
    -------
    pd.DataFrame or None
        blocks.det DataFrame, or None if no blocks found.
    """
    if outdir is None:
        outdir = os.path.join(os.path.dirname(bfile), "_ld_tmp")
    os.makedirs(outdir, exist_ok=True)

    keep_args = []
    if sample_file:
        keep_args = ["--keep", sample_file]

    prefix = os.path.join(outdir, f"blocks_chr{chr_num}")
    cmd = [plink_path, "--bfile", bfile] + keep_args + [
        "--chr", str(chr_num), "--blocks", "no-pheno-req",
        "--blocks-max-kb", str(blocks_max_kb),
        "--blocks-min-maf", str(blocks_min_maf),
        "--out", prefix,
    ]
    subprocess.run(cmd, capture_output=True, text=True)

    det_file = f"{prefix}.blocks.det"
    if os.path.exists(det_file) and os.path.getsize(det_file) > 0:
        return pd.read_csv(det_file, sep=r'\s+')
    return None


def compute_ld_decay(
    bfile,
    chr_num,
    plink_path="plink",
    sample_file=None,
    ld_window=200,
    ld_window_kb=1000,
    bin_size_bp=10000,
    max_distance_bp=500000,
    outdir=None,
):
    """
    Compute LD decay curve for one chromosome.

    Parameters
    ----------
    bfile : str
        PLINK binary file prefix.
    chr_num : int
        Chromosome number.
    bin_size_bp : int
        Distance bin size for averaging r².
    max_distance_bp : int
        Maximum distance to consider.

    Returns
    -------
    pd.DataFrame
        Columns: ['distance_bp', 'mean_r2', 'n_pairs']
    """
    if outdir is None:
        outdir = os.path.join(os.path.dirname(bfile), "_ld_tmp")
    os.makedirs(outdir, exist_ok=True)

    keep_args = []
    if sample_file:
        keep_args = ["--keep", sample_file]

    prefix = os.path.join(outdir, f"ld_chr{chr_num}")
    cmd = [plink_path, "--bfile", bfile] + keep_args + [
        "--chr", str(chr_num),
        "--r2", "--ld-window", str(ld_window),
        "--ld-window-kb", str(ld_window_kb),
        "--ld-window-r2", "0",
        "--out", prefix,
    ]
    subprocess.run(cmd, capture_output=True, text=True)

    ld_file = f"{prefix}.ld"
    if not os.path.exists(ld_file):
        return None

    # Parse and bin
    bins = np.arange(0, max_distance_bp + bin_size_bp, bin_size_bp)
    sums = np.zeros(len(bins) - 1)
    counts = np.zeros(len(bins) - 1, dtype=int)

    with open(ld_file) as f:
        header = f.readline()
        for line in f:
            parts = line.split()
            if len(parts) >= 7:
                bp_a = int(parts[1])
                bp_b = int(parts[4])
                r2 = float(parts[6])
                dist = abs(bp_b - bp_a)
                idx = int(dist // bin_size_bp)
                if 0 <= idx < len(sums):
                    sums[idx] += r2
                    counts[idx] += 1

    mask = counts > 0
    result = pd.DataFrame({
        "distance_bp": (bins[:-1][mask] + bin_size_bp / 2).astype(int),
        "mean_r2": sums[mask] / counts[mask],
        "n_pairs": counts[mask],
    })
    return result
