"""
BFI core computation engine.

    BFI_i = sum_b RAB_ib * w_b / sum_b w_b
    w_b   = max(r2_obs,b - r2_eq,b, 0)                      (excess-LD weight)

r2_obs,b  mean squared Pearson correlation of additive genotype codes over the SNP pairs P_b of
          block b. r2_mode="pairwise" (default, v1.0.2): each pair uses the reference accessions
          genotyped at both SNPs (pairwise deletion, as PLINK --r2), pairs observed in < 50% of the
          reference accessions or monomorphic in their overlap are skipped. r2_mode="complete"
          (v1.0.0-1.0.1): only SNPs with no missing genotype in the whole reference group are used,
          and blocks with < 2 such SNPs are dropped.
r2_eq,b   drift-recombination expectation 1/(1 + 4 Ne c d) (Sved 1971; Hill 1981):
            eq_mode="span" (default, the published definition): d = block span BP2 - BP1 from
                           PLINK --blocks, i.e. the equilibrium LD expected across the whole block
            eq_mode="pair" (sensitivity option): mean over the SNP pairs P_b of the
                           expectation at each pair distance |x_i - x_j|
          optionally + 1/n (sampling term).
RAB_ib    proportion of non-missing block SNPs at which accession i is not homozygous for the
          reference-group major allele (heterozygotes count as non-conforming).

With `reference_file`, LD blocks, major alleles and weights come from the reference samples and
every sample is scored against them; all samples are extracted ONCE so that SNP identity, order
and allele coding are shared by reference and scored samples.
"""

import os
import subprocess
import shutil
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

from bfi.ne import equilibrium_r2

EQ_MODES = ("span", "pair")
R2_MODES = ("pairwise", "complete")


# ---------------------------------------------------------------- pure functions (unit-tested)
def block_r2_obs(mat, r2_mode="pairwise", min_frac=0.5):
    """Mean squared correlation over the SNP pairs of a block (rows = accessions, NaN = missing).

    Returns (r2_obs, pair_i, pair_j) with the column indices of the averaged SNP pairs,
    or (None, None, None) if no pair is usable."""
    k = mat.shape[1]
    if r2_mode == "complete":
        is_valid = ~np.isnan(mat)
        use = np.where(np.all(is_valid, axis=0))[0]
        if len(use) < 2:
            return None, None, None
        sub = mat[:, use]
        sd = np.std(sub, axis=0)
        good = sd > 0
        if good.sum() < 2:
            return None, None, None
        z = (sub[:, good] - sub[:, good].mean(axis=0)) / sd[good]
        corr = np.corrcoef(z.T)
        a, b = np.triu_indices(corr.shape[0], k=1)
        cols = use[good]
        return float(np.mean(corr[a, b] ** 2)), cols[a], cols[b]
    if r2_mode != "pairwise":
        raise ValueError(f"r2_mode must be one of {R2_MODES}")
    M = ~np.isnan(mat)
    X = np.where(M, mat, 0.0)
    Mf = M.astype(np.float64)
    nij = Mf.T @ Mf
    S1 = X.T @ Mf
    S2 = (X * X).T @ Mf
    Sxy = X.T @ X
    num = nij * Sxy - S1 * S1.T
    vi = nij * S2 - S1 ** 2
    a, b = np.triu_indices(k, k=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = num[a, b] / np.sqrt(vi[a, b] * vi.T[a, b])
    ok = (nij[a, b] >= min_frac * mat.shape[0]) & np.isfinite(r) & (vi[a, b] > 0) & (vi.T[a, b] > 0)
    if not ok.any():
        return None, None, None
    return float(np.mean(r[ok] ** 2)), a[ok], b[ok]


def block_r2_eq(span_bp, pair_dist, ne, recomb_rate, eq_mode="span", n=None):
    """Expected r2 of a block.

    span_bp   : block span BP2 - BP1 as reported by PLINK --blocks (used by eq_mode="span")
    pair_dist : distances (bp) of the SNP pairs whose r2 was averaged (used by eq_mode="pair")"""
    if eq_mode == "span":
        e = float(equilibrium_r2(span_bp, ne, recomb_rate))
    elif eq_mode == "pair":
        e = float(np.mean(equilibrium_r2(np.asarray(pair_dist, dtype=np.float64), ne, recomb_rate)))
    else:
        raise ValueError(f"eq_mode must be one of {EQ_MODES}")
    if n:
        e += 1.0 / n
    return e


def major_homozygote(mat_ref):
    """Dosage code (0 or 2) of the reference-group major homozygote per SNP."""
    return np.where(np.nanmean(mat_ref, axis=0) <= 1.0, 0, 2)


def rare_allele_burden(mat, major):
    """Per-row proportion of non-missing SNPs that are not homozygous for `major`."""
    is_valid = ~np.isnan(mat)
    matches = (mat == major[np.newaxis, :]) & is_valid
    vc = is_valid.sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(vc > 0, 1.0 - matches.sum(axis=1) / np.maximum(vc, 1), 0.0)


# ---------------------------------------------------------------- driver
def check_unique_ids(bfile):
    ids = pd.read_csv(f"{bfile}.bim", sep=r"\s+", header=None, usecols=[1], dtype=str)[1]
    dup = ids.duplicated().sum()
    if dup:
        raise ValueError(f"{bfile}.bim contains {dup} duplicated variant IDs; remove duplicates first "
                         f"(e.g. plink2 --rm-dup exclude-mismatch --make-bed).")


def compute_bfi(
    bfile,
    sample_file=None,
    ne=None,
    recomb_rate=4e-8,
    chromosomes=None,
    blocks_max_kb=200,
    blocks_min_maf=0.05,
    maf=None,
    min_block_snps=2,
    eq_mode="span",
    r2_mode="pairwise",
    sampling_term=False,
    reference_file=None,
    n_workers=6,
    plink_path="plink",
    tmpdir=None,
    return_blocks=False,
):
    """
    Compute genome-wide BFI.

    Parameters
    ----------
    bfile : str             PLINK binary prefix (.bed/.bim/.fam); variant IDs must be unique.
    sample_file : str       Samples to score (IID per line, or FID<TAB>IID). Default: all.
    ne : float              Effective population size of the reference group (required).
    recomb_rate : float     Per-bp recombination rate (rice 4e-8, soybean 2.5e-8, maize 1e-8).
    blocks_max_kb, blocks_min_maf : PLINK --blocks settings.
    maf : float, optional   Extra --maf filter applied with --blocks (maize analyses: 0.05).
    eq_mode : "span"|"pair" Expected-r2 definition (see module docstring); "span" = published.
    r2_mode : "pairwise"|"complete" Missing-data handling for r2_obs (see module docstring).
    sampling_term : bool    Add 1/n to r2_eq.
    reference_file : str    Reference samples defining blocks, major alleles and weights
                            (cross-group BFI). Default: the scored samples themselves.
    return_blocks : bool    Also return the block table.

    Returns
    -------
    DataFrame ['IID', 'BFI'] (sorted by BFI, descending) [, block table]
    """
    if ne is None:
        raise ValueError("ne (effective population size) is required")
    if eq_mode not in EQ_MODES:
        raise ValueError(f"eq_mode must be one of {EQ_MODES}")
    if r2_mode not in R2_MODES:
        raise ValueError(f"r2_mode must be one of {R2_MODES}")
    check_unique_ids(bfile)
    if chromosomes is None:
        bim = pd.read_csv(f"{bfile}.bim", sep=r"\s+", header=None, usecols=[0], dtype=str)
        chromosomes = list(dict.fromkeys(bim[0]))
    if tmpdir is None:
        tmpdir = os.path.join(os.path.dirname(os.path.abspath(bfile)), "_bfi_tmp")
    args = [(bfile, sample_file, reference_file, ne, recomb_rate, c, blocks_max_kb, blocks_min_maf,
             maf, min_block_snps, eq_mode, sampling_term, plink_path, tmpdir, r2_mode) for c in chromosomes]
    res = []
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        for fut in as_completed([ex.submit(compute_bfi_chromosome, *a) for a in args]):
            r = fut.result()
            if r is not None:
                res.append(r)
    if not res:
        raise RuntimeError("No blocks found on any chromosome")
    base = res[0][0][["IID"]].copy()
    for r in res:
        assert np.array_equal(r[0]["IID"].values, base["IID"].values)
    wb = sum(r[0]["w_bfi"].values for r in res)
    wt = sum(r[0]["w_total"].values[0] for r in res)
    base["BFI"] = wb / wt if wt > 0 else 0.0
    out = base.sort_values("BFI", ascending=False).reset_index(drop=True)
    if return_blocks:
        return out, pd.concat([r[1] for r in res], ignore_index=True)
    return out


def _keep(path, dst):
    with open(path) as fin, open(dst, "w") as f:
        for line in fin:
            p = line.split()
            if p:
                f.write(f"{p[0]}\t{p[1]}\n" if len(p) > 1 else f"0\t{p[0]}\n")


def compute_bfi_chromosome(bfile, sample_file, reference_file, ne, recomb_rate, chr_num,
                           blocks_max_kb, blocks_min_maf, maf, min_block_snps, eq_mode,
                           sampling_term, plink_path, tmpdir_base, r2_mode="pairwise"):
    tmpdir = os.path.join(tmpdir_base, f"chr{chr_num}_{os.getpid()}")
    os.makedirs(tmpdir, exist_ok=True)
    try:
        return _compute_chr_inner(bfile, sample_file, reference_file, ne, recomb_rate, chr_num,
                                  blocks_max_kb, blocks_min_maf, maf, min_block_snps, eq_mode,
                                  sampling_term, plink_path, tmpdir, r2_mode)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def _compute_chr_inner(bfile, sample_file, reference_file, ne, recomb_rate, chr_num,
                       blocks_max_kb, blocks_min_maf, maf, min_block_snps, eq_mode,
                       sampling_term, plink_path, tmpdir, r2_mode="pairwise"):
    ref = reference_file or sample_file
    ref_args = []
    if ref:
        _keep(ref, os.path.join(tmpdir, "ref.txt"))
        ref_args = ["--keep", os.path.join(tmpdir, "ref.txt")]
    # 1. LD blocks within the reference samples
    bp = os.path.join(tmpdir, "blocks")
    cmd = [plink_path, "--bfile", bfile] + ref_args + [
        "--chr", str(chr_num), "--blocks", "no-pheno-req",
        "--blocks-max-kb", str(blocks_max_kb), "--blocks-min-maf", str(blocks_min_maf)]
    if maf is not None:
        cmd += ["--maf", str(maf)]
    subprocess.run(cmd + ["--out", bp], capture_output=True, text=True)
    bf = f"{bp}.blocks.det"
    if not os.path.exists(bf) or os.path.getsize(bf) == 0:
        return None
    blocks = pd.read_csv(bf, sep=r"\s+")
    if len(blocks) == 0:
        return None
    snps = sorted({s for x in blocks["SNPS"] for s in x.split("|")})
    sf = os.path.join(tmpdir, "snps.txt")
    with open(sf, "w") as f:
        f.write("\n".join(snps) + "\n")
    # 2. ONE extraction of reference + scored samples (shared SNP order and allele coding)
    score_args = []
    if sample_file and reference_file:
        _keep(sample_file, os.path.join(tmpdir, "s.txt"))
        both = pd.concat([pd.read_csv(os.path.join(tmpdir, f), sep="\t", header=None, dtype=str)
                          for f in ("ref.txt", "s.txt")]).drop_duplicates()
        both.to_csv(os.path.join(tmpdir, "both.txt"), sep="\t", header=False, index=False)
        score_args = ["--keep", os.path.join(tmpdir, "both.txt")]
    elif sample_file:
        score_args = ref_args
    elif reference_file:
        score_args = []           # score everybody
    gp = os.path.join(tmpdir, "geno")
    subprocess.run([plink_path, "--bfile", bfile] + score_args +
                   ["--extract", sf, "--chr", str(chr_num), "--recode", "A", "--out", gp],
                   capture_output=True, text=True)
    if not os.path.exists(f"{gp}.raw"):
        return None
    geno = pd.read_csv(f"{gp}.raw", sep=r"\s+")
    cols = [c for c in geno.columns if c not in ("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")]
    snp_of = {"_".join(c.split("_")[:-1]): c for c in cols}
    bim = pd.read_csv(f"{bfile}.bim", sep=r"\s+", header=None, usecols=[0, 1, 3], names=["chr", "id", "bp"],
                      dtype={"chr": str, "id": str, "bp": np.int64})
    bim = bim[bim.chr == str(chr_num)]
    pos_of = dict(zip(bim.id, bim.bp))
    is_ref = np.ones(len(geno), dtype=bool)
    if ref:
        ref_ids = set(pd.read_csv(os.path.join(tmpdir, "ref.txt"), sep="\t", header=None, dtype=str)[1])
        is_ref = geno["IID"].astype(str).isin(ref_ids).values
    score_rows = np.ones(len(geno), dtype=bool)
    if sample_file:
        s_ids = set(pd.read_csv(os.path.join(tmpdir, "s.txt" if reference_file else "ref.txt"),
                                sep="\t", header=None, dtype=str)[1])
        score_rows = geno["IID"].astype(str).isin(s_ids).values
    n_ref = int(is_ref.sum())
    G = geno[cols].values.astype(np.float64)
    Gref, Gs = G[is_ref], G[score_rows]
    col_idx = {c: k for k, c in enumerate(cols)}
    w_bfi = np.zeros(int(score_rows.sum()))
    w_total, rows = 0.0, []
    for _, row in blocks.iterrows():
        bc = [snp_of[s] for s in row["SNPS"].split("|") if s in snp_of]
        if len(bc) < min_block_snps:
            continue
        idx = [col_idx[c] for c in bc]
        mref = Gref[:, idx]
        r2o, pi, pj = block_r2_obs(mref, r2_mode)
        if r2o is None:
            continue
        pos = np.array([pos_of["_".join(c.split("_")[:-1])] for c in bc], dtype=np.float64)
        r2e = block_r2_eq(int(row["BP2"]) - int(row["BP1"]), np.abs(pos[pj] - pos[pi]), ne, recomb_rate,
                          eq_mode, n_ref if sampling_term else None)
        w = max(r2o - r2e, 0.0)
        rows.append((chr_num, int(row["BP1"]), int(row["BP2"]), len(bc), len(pi), r2o, r2e, w))
        if w <= 0:
            continue
        rab = rare_allele_burden(Gs[:, idx], major_homozygote(mref))
        w_bfi += rab * w
        w_total += w
    out = pd.DataFrame({"IID": geno["IID"].values[score_rows], "w_bfi": w_bfi, "w_total": w_total})
    btab = pd.DataFrame(rows, columns=["chr", "bp1", "bp2", "n_snps", "n_pairs_r2", "r2_obs", "r2_eq", "w"])
    return out, btab
