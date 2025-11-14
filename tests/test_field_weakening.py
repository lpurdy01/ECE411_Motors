import math

import numpy as np

from motor_models import (
    FieldWeakeningParams,
    constant_power_torque,
    field_weakening_characteristics,
    field_weakening_flux,
)


def test_field_weakening_flux_matches_hand_calculation():
    v0 = 1.0
    i0 = 0.5
    ls = 0.3
    w_e = 1.2
    expected = math.sqrt((v0 / w_e) ** 2 + (ls * i0) ** 2)
    assert math.isclose(field_weakening_flux(v0, i0, ls, w_e), expected, rel_tol=1e-9)


def test_constant_power_relationships():
    te, p = constant_power_torque(1.1, 0.9, 1.4)
    assert math.isclose(te * 1.4, p, rel_tol=1e-9)


def test_characteristics_monotonic_flux():
    params = FieldWeakeningParams(v0_max=1.0, i0_max=1.2, ls=0.5, rs=0.0)
    curves = field_weakening_characteristics(params, np.linspace(0.4, 2.0, 20))
    assert np.all(np.diff(curves["flux"]) <= 1e-9)
