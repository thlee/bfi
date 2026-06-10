# BFI — Breeding Freedom Index

A quantitative metric for genomic design space in crop breeding.

BFI measures an individual accession's potential to disrupt excess LD blocks when used as a crossing parent. It weights each LD block by its excess LD above the drift–recombination equilibrium baseline (Hill, 1981), enabling cross-species comparison without requiring wild-type reference samples.

## Installation

```bash
pip install bfi
```

**Requirement:** [PLINK 1.9](https://www.cog-genomics.org/plink/) must be installed and accessible.

## Quick Start

### Command line

```bash
# Compute BFI for rice indica (N_e = 911, recombination rate = 4 cM/Mb)
bfi run --bfile my_data --ne 911 --recomb-rate 4e-8 --out bfi_results.tsv

# Specify chromosomes and number of workers
bfi run --bfile my_data --ne 911 --chromosomes 1,2,3 --workers 12 --out bfi_chr1-3.tsv

# Custom PLINK path
bfi run --bfile my_data --ne 911 --plink /path/to/plink --out results.tsv
```

### Python API

```python
from bfi import compute_bfi

result = compute_bfi(
    bfile="my_data",          # PLINK .bed/.bim/.fam prefix
    ne=911,                    # Effective population size
    recomb_rate=4e-8,          # Per-bp recombination rate
    chromosomes=[1, 2, 3],     # Optional: specific chromosomes
    n_workers=6,               # Parallel workers
    plink_path="plink",        # PLINK 1.9 path
)

# result is a DataFrame with columns ['IID', 'BFI']
print(result.head())
```

## Parameters

| Parameter | Description | Typical values |
|-----------|-------------|----------------|
| `ne` | Effective population size | Rice indica: 911, Soybean improved: 813 |
| `recomb_rate` | Per-bp recombination rate | Rice: 4e-8, Soybean: 2.5e-8, Maize: 1e-8 |
| `blocks_max_kb` | Max LD block size (kb) | 200 (robust to 50–500) |
| `blocks_min_maf` | Min MAF for block SNPs | 0.05 |

## How BFI Works

```
BFI = Σ_b [ RAB_b × w_b ] / Σ_b [ w_b ]

RAB_b = rare allele burden in block b
w_b   = max(r²_obs - r²_eq, 0)   ← excess LD weight
r²_eq = 1 / (1 + 4·N_e·c)        ← Hill (1981) equilibrium
```

- **BFI ≈ 0**: Conforms to population's LD-locked haplotypes → rigid, low breeding freedom
- **BFI ≈ 1**: Divergent at high-LD blocks → high potential for LD disruption (Break Shot candidate)

## Citation

```
Lee, T.-H. (2026). Breeding Freedom Index: A Quantitative Metric for
Genomic Design Space in Crop Breeding.
```

## License

MIT
