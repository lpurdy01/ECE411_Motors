import math

import numpy as np

from motor_models import (
    abc_to_alphabeta,
    abc_to_qd0,
    alphabeta_to_abc,
    alphabeta_to_qd,
    qd_to_alphabeta,
)


def test_balanced_three_phase_transforms_align_with_q_axis():
    theta = 0.0
    va = np.cos(theta)
    vb = np.cos(theta - 2 * math.pi / 3)
    vc = np.cos(theta - 4 * math.pi / 3)

    alpha, beta, v0 = abc_to_alphabeta(va, vb, vc)
    assert abs(alpha - 1.0) < 1e-12
    assert abs(beta) < 1e-12
    assert abs(v0) < 1e-12

    vq, vd, v0_again = abc_to_qd0(np.array([va]), np.array([vb]), np.array([vc]), np.array([0.0]))
    assert np.allclose(vq, 1.0)
    assert np.allclose(vd, 0.0)
    assert np.allclose(v0_again, 0.0)


def test_quadrature_rotation_swaps_q_and_d():
    alpha = np.array([1.0])
    beta = np.array([0.0])
    theta = np.array([math.pi / 2])

    vq, vd = alphabeta_to_qd(alpha, beta, theta)
    assert np.allclose(vq, 0.0, atol=1e-12)
    assert np.allclose(vd, 1.0, atol=1e-12)


def test_round_trip_qd_to_abc_returns_original_phases():
    rng = np.random.default_rng(1234)
    vq = rng.normal(size=400)
    vd = rng.normal(size=400)
    theta = np.linspace(0.0, 2 * math.pi, 400)
    alpha, beta = qd_to_alphabeta(vq, vd, theta)
    va, vb, vc = alphabeta_to_abc(alpha, beta, np.zeros_like(alpha))
    recovered_vq, recovered_vd, _ = abc_to_qd0(va, vb, vc, theta)
    assert np.allclose(recovered_vq, vq, atol=1e-9)
    assert np.allclose(recovered_vd, vd, atol=1e-9)
