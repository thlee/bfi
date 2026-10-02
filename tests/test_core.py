"""Unit tests for the BFI building blocks (no PLINK needed): python -m pytest tests/"""
import numpy as np
import pandas as pd
import pytest

from bfi.core import block_r2_obs, block_r2_eq, major_homozygote, rare_allele_burden
from bfi.ne import equilibrium_r2, estimate_ne


def test_equilibrium_r2():
    assert equilibrium_r2(0, 1000, 4e-8) == 1.0
    assert np.isclose(equilibrium_r2(10_000, 1000, 4e-8), 1 / (1 + 4 * 1000 * 4e-4))


def test_span_and_pair_expectations():
    pos = np.array([0, 100, 250, 5_000, 40_000], dtype=float)
    i, j = np.triu_indices(len(pos), 1)
    d = np.abs(pos[j] - pos[i])
    span = block_r2_eq(40_000, d, 911, 4e-8, "span")
    pair = block_r2_eq(None, d, 911, 4e-8, "pair")
    assert np.isclose(span, 1 / (1 + 4 * 911 * 4e-8 * 40_000))      # block span from PLINK BP1/BP2
    assert np.isclose(pair, np.mean(1 / (1 + 4 * 911 * 4e-8 * d)))
    assert pair > span                                               # E[r2] decreases with distance
    assert np.isclose(block_r2_eq(None, d, 911, 4e-8, "pair", n=100), pair + 0.01)


def test_r2_obs_complete_mode_uses_complete_polymorphic_snps_only():
    g = np.array([[0, 0, 2, 0], [2, 2, 0, 0], [0, 0, 2, np.nan], [2, 2, 0, 0]], dtype=float)
    r2, pi, pj = block_r2_obs(g, "complete")
    assert set(pi) | set(pj) == {0, 1, 2}                            # SNP4 has a missing call
    assert np.isclose(r2, 1.0)


def test_r2_obs_pairwise_matches_pandas():
    rng = np.random.default_rng(7)
    n, k = 200, 8
    base = rng.integers(0, 3, size=(n, 1)).astype(float)
    m = np.where(rng.random((n, k)) < 0.7, base, rng.integers(0, 3, size=(n, k))).astype(float)
    m[rng.random((n, k)) < 0.05] = np.nan
    r2, pi, pj = block_r2_obs(m, "pairwise")
    c = pd.DataFrame(m).corr(min_periods=100).values[np.triu_indices(k, 1)]
    assert np.isclose(r2, np.nanmean(c ** 2))
    assert len(pi) == np.isfinite(c).sum()


def test_pairwise_keeps_blocks_that_complete_mode_drops():
    g = np.array([[0, 0], [2, 2], [0, np.nan], [2, 2], [np.nan, 0], [0, 0]], dtype=float)
    assert block_r2_obs(g, "complete")[0] is None                    # no SNP is complete
    r2, _, _ = block_r2_obs(g, "pairwise")
    assert np.isclose(r2, 1.0)


def test_rab_counts_heterozygotes_as_nonconforming():
    ref = np.array([[0, 0], [0, 0], [0, 2], [2, 0]], dtype=float)
    major = major_homozygote(ref)
    assert major.tolist() == [0, 0]
    x = np.array([[0, 0], [1, 0], [2, 2], [np.nan, 1]], dtype=float)
    assert np.allclose(rare_allele_burden(x, major), [0.0, 0.5, 1.0, 1.0])


def test_estimate_ne_recovers_truth():
    d = np.arange(5_000, 1_000_000, 10_000, dtype=float)
    r2 = 1 / (1 + 4 * 900 * 4e-8 * d) + 1 / 200
    fit = estimate_ne(d, r2, recomb_rate=4e-8, n_samples=200)
    assert abs(fit["ne"] - 900) < 1.0


def test_bad_modes():
    with pytest.raises(ValueError):
        block_r2_eq(10, np.array([10.0]), 100, 1e-8, "midpoint")
    with pytest.raises(ValueError):
        block_r2_obs(np.zeros((3, 2)), "listwise")
