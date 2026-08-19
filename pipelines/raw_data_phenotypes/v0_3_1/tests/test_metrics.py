import numpy as np

from neurothermo_phenotypes.metrics import (
    block_path_kl_rate,
    electrical_work,
    permutation_entropy,
)


def test_external_work_units():
    t = np.linspace(0, 1, 1001)
    v = np.full_like(t, -60.0)
    result = electrical_work(t, v, 100.0, -70.0)
    assert abs(result["external_work_signed_fJ"] - 1000.0) < 1e-6


def test_permutation_entropy_orders_noise_above_constant():
    rng = np.random.default_rng(2)
    constant = np.ones(1000)
    noise = rng.normal(size=1000)
    assert permutation_entropy(constant, 4, 1) == 0.0
    assert permutation_entropy(noise, 4, 1) > 0.9


def test_path_kl_detects_asymmetric_cycle():
    phase = np.linspace(0, 50, 8000)
    asymmetric = (phase % 1.0) ** 4
    reversible = np.sin(2 * np.pi * phase)
    d_asym, _, _ = block_path_kl_rate(asymmetric, 0.001, n_bins=6, word_length=3, delay=2)
    d_rev, _, _ = block_path_kl_rate(reversible, 0.001, n_bins=6, word_length=3, delay=2)
    assert d_asym > d_rev

