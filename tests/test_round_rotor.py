import math

import numpy as np

from motor_models import RoundRotorParams, round_rotor_operating_point


def test_round_rotor_solution_satisfies_voltage_equations():
    params = RoundRotorParams(rs=0.05, ls=0.7, v0=1.05, e0=0.95, w_e=1.0)
    delta = math.radians(25)
    op = round_rotor_operating_point(params, delta)

    vq = params.v0 * math.cos(delta)
    vd = params.v0 * math.sin(delta)
    lhs_q = params.rs * op["iq"] + params.w_e * params.ls * op["id"] + op["eq"]
    lhs_d = params.rs * op["id"] - params.w_e * params.ls * op["iq"] + op["ed"]

    assert math.isclose(lhs_q, vq, rel_tol=1e-9, abs_tol=1e-9)
    assert math.isclose(lhs_d, vd, rel_tol=1e-9, abs_tol=1e-9)

    expected_power = 1.5 * (vq * op["iq"] + vd * op["id"])
    assert math.isclose(expected_power, op["power"], rel_tol=1e-9)


def test_round_rotor_zero_delta_yields_quadrature_current_zero():
    params = RoundRotorParams(rs=0.0, ls=0.6, v0=1.0, e0=0.8, w_e=1.0)
    op = round_rotor_operating_point(params, delta=0.0)
    assert math.isclose(op["iq"], 0.0, abs_tol=1e-12)
    assert op["power"] == 0.0
