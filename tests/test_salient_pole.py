import math

import numpy as np

from motor_models import SalientPoleParams, salient_pole_torque_curves


def test_salient_torque_matches_closed_form_extrema():
    params = SalientPoleParams(ld=1.2, lq=0.8, e0=1.0, v0=1.0, w_e=1.0)
    delta = np.linspace(-math.pi, math.pi, 600)
    t_field, t_rel, total = salient_pole_torque_curves(params, delta)

    expected_field_max = (params.e0 * params.v0) / (params.w_e * params.ld)
    expected_rel_max = abs(params.v0**2 * (params.ld - params.lq) / (2 * params.w_e * params.ld * params.lq))

    assert math.isclose(np.max(np.abs(t_field)), expected_field_max, rel_tol=1e-3)
    assert math.isclose(np.max(np.abs(t_rel)), expected_rel_max, rel_tol=1e-3)
    assert np.allclose(total, t_field + t_rel)


def test_salient_symmetry_about_zero_delta():
    params = SalientPoleParams(ld=1.1, lq=0.9, e0=0.95, v0=1.02, w_e=1.0)
    delta = np.array([0.0, math.pi / 4, -math.pi / 4])
    t_field, t_rel, total = salient_pole_torque_curves(params, delta)
    assert math.isclose(t_field[1], -t_field[2], rel_tol=1e-12, abs_tol=1e-12)
    assert math.isclose(t_rel[1], -t_rel[2], rel_tol=1e-12, abs_tol=1e-12)
    assert math.isclose(total[1], -total[2], rel_tol=1e-12, abs_tol=1e-12)
    assert total[0] == 0.0
