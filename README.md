# BFI — Breeding Freedom Index

BFI ranks accessions by their rare-allele burden in LD blocks, weighting each block by its
excess LD above the drift–recombination expectation (Hill 1981). It needs only the genotype
data of a reference group and that group's effective population size (Nₑ).

## Installation

```bash
pip install git+https://github.com/thlee/bfi.git@v1.0.2
# archived release: Zenodo, DOI: TBD
```

The name `bfi` on PyPI belongs to an unrelated project — do **not** use `pip install bfi`.
**Requirement:** [PLINK 1.9](https://www.cog-genomics.org/plink/) on the PATH (or `--plink`).
Variant IDs in the `.bim` must be unique (the tool stops otherwise; remove duplicates with
`plink2 --rm-dup exclude-mismatch --make-bed`).

## Quick start

```bash
# within-group BFI (rice indica, Ne = 911, 4 cM/Mb)
bfi run --bfile rice --keep indica.txt --ne 911 --recomb-rate 4e-8 --out bfi_indica.tsv

# cross-group BFI: every accession scored against the improved soybean reference
bfi run --bfile soybean --reference improved.txt --ne 1265 --recomb-rate 2.5e-8 --out bfi_cross.tsv

# maize settings used in the article (block SNPs MAF >= 0.10 after a 0.05 MAF filter)
bfi run --bfile maize --keep improved.txt --ne 42846 --recomb-rate 1e-8 \
        --blocks-min-maf 0.10 --maf 0.05 --blocks-out blocks.tsv --out bfi_maize.tsv
```

```python
from bfi import compute_bfi
res = compute_bfi(bfile="rice", sample_file="indica.txt", ne=911, recomb_rate=4e-8)
```

## Definition

```
BFI_i   = Σ_b RAB_ib · w_b / Σ_b w_b
w_b     = max(r²_obs,b − r²_eq,b, 0)
r²_eq,b = 1 / (1 + 4·Nₑ·c·(BP2 − BP1))                                       (--eq-mode span, default)
          mean over the SNP pairs of r²_obs,b of 1 / (1 + 4·Nₑ·c·d_ij)       (--eq-mode pair, sensitivity)
          [+ 1/n with --sampling-term]
```

The default compares the block's mean r² with the equilibrium r² expected across the whole block
span (Gabriel blocks keep strong LD up to their ends), so the weight measures LD that persists over
the block's physical extent. `--eq-mode pair` evaluates the expectation at every SNP-pair distance
instead; with the simple Sved/Hill expectation this gives nearly identical accession rankings in the
article's data but down-weights short blocks further.

* Blocks: PLINK `--blocks` (Gabriel et al. 2002) within the reference samples.
* r²_obs,b: mean squared correlation of additive genotypes over the block's SNP pairs. Default
  `--r2-mode pairwise`: each pair uses the reference accessions genotyped at both SNPs (as PLINK
  `--r2`; pairs observed in < 50% of the reference accessions are skipped). `--r2-mode complete`
  reproduces v1.0.0–1.0.1 (only SNPs without any missing call; blocks with < 2 such SNPs dropped),
  which discards most blocks in large, heterogeneously sequenced panels.
* RAB_ib: share of accession *i*'s non-missing block SNPs that are not homozygous for the
  reference-group major allele (heterozygotes count as non-conforming).
* With `--reference`, the reference and scored samples are extracted together, so SNP identity,
  order and allele coding are shared (fixes the A1-recoding issue of separate PLINK extractions).

| Parameter | Meaning | Values used in the article |
|---|---|---|
| `--ne` | Nₑ of the reference group (Hill fit to its LD decay) | Table 1 |
| `--recomb-rate` | per-bp recombination rate | rice 4e-8, soybean 2.5e-8, maize 1e-8 |
| `--blocks-max-kb` | max block length | 200 |
| `--blocks-min-maf` / `--maf` | block SNP filters | 0.05 / – (rice, soybean); 0.10 / 0.05 (maize) |

`bfi.estimate_ne(distances, r2, recomb_rate, n_samples)` fits Nₑ to binned LD decay
(r² − 1/n correction when `n_samples` is given).

## Changes

* **1.0.2** — r²_obs with pairwise-complete data (`--r2-mode pairwise`, new default; `complete`
  keeps the 1.0.1 behaviour); `--reference` for cross-group scoring with one shared extraction and
  explicit allele alignment (fixes SNPs silently dropped when PLINK re-assigned A1 between separate
  extractions); `--eq-mode pair` and `--sampling-term` sensitivity options (default `span` unchanged);
  `--maf`; duplicate-variant-ID check; block table output; Nₑ fit options aligned with the analyses;
  unit tests.
* 1.0.1 — licence/citation for five authors, environment files.
* 1.0.0 — first public release.

## Citation

Kim M, Lee K, Kim KD, Hwang J-H, Lee T-H (2026) Breeding Freedom Index: an
excess-linkage-disequilibrium-weighted framework for genomic design space in crop breeding.
(Journal details to be added on publication.) Software/data archive: Zenodo, DOI: TBD.

## License

MIT
