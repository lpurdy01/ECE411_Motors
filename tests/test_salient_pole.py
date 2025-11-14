import math

import numpy as np

from motor_models import SalientPoleParams, salient_pole_torque_curves


def test_salient_torque_matches_closed_form_extrema():
    params = SalientPoleParams(ld=1.2, lq=0.8, e0=1.0, v0=1.0, w_e=1.0)
    delta = np.linspace(-math.pi, math.pi, 600)
    t_field, t_rel, _ = salient_pole_torque_curves(params, delta)

    expected_field_max = (params.e0 * params.v0) / (params.w_e * params.ld)
    expected_rel_max = abs(params.v0**2 * (params.ld - params.lq) / (2 * params.w_e * params.ld * params.lq))

    assert math.isclose(np.max(np.abs(t_field)), expected_field_max, rel_tol=1e-3)
    assert math.isclose(np.max(np.abs(t_rel)), expected_rel_max, rel_tol=1e-3)
