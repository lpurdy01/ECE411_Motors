import math

import numpy as np

from motor_models import (
    default_machine_params,
    inverter_voltage_vs_speed,
    mechanical_power_map,
    pu_bases_from_machine,
    sm_required_voltage_park,
    sm_required_voltage_park_pu,
)


def test_sm_required_voltage_matches_manual():
    params = default_machine_params()
    omega_m = 200.0  # rad/s mechanical
    iq = 150.0
    id_ = 0.0
    field_current = 5.0

    result = sm_required_voltage_park(omega_m, field_current, iq, id_, params)

    w_e = (params.p / 2.0) * omega_m
    expected_v = (w_e * params.Lm * field_current) + (params.Rs + 1j * w_e * params.Ls) * (iq + 0j)
    assert math.isclose(result["w_e"], w_e)
    assert np.isclose(result["V"], expected_v).all()
    assert math.isclose(result["V0"], abs(expected_v))


def test_pu_wrapper_and_bases_consistent():
    params = default_machine_params()
    bases = pu_bases_from_machine(params)

    omega_m = np.array([0.0, params.wmb])
    iq = np.array([0.0, params.I0_max])
    res_pu = sm_required_voltage_park_pu(omega_m, params.If_max, iq, 0.0, params, bases)
    res_si = sm_required_voltage_park(omega_m, params.If_max, iq, 0.0, params)

    # Base voltage follows Vdc/sqrt(3)
    assert math.isclose(bases.VB, params.Vdc / math.sqrt(3.0), rel_tol=1e-6)
    # Per-unit magnitudes convert back to SI via VB
    assert np.allclose(res_si["V0"], res_pu["V0_pu"] * bases.VB)


def test_inverter_voltage_crossing_detected():
    params = default_machine_params()
    bases = pu_bases_from_machine(params)
    speeds = np.linspace(0.0, 3 * params.wmb, 50)
    sweep = inverter_voltage_vs_speed(params, bases, [params.If_max], speeds)

    trace = sweep["traces"][0]
    assert trace["cross_speed"] is not None
    assert trace["V0"].shape == speeds.shape


def test_power_map_masks_infeasible_points():
    params = default_machine_params()
    bases = pu_bases_from_machine(params)
    torque_axis = np.linspace(0.0, 250.0, 20)
    omega_axis = np.linspace(0.0, 200.0, 20)

    power_map = mechanical_power_map(params.If_max, torque_axis, omega_axis, params, bases)
    # At zero speed (first column), any nonzero torque should be infeasible
    assert not power_map["feasible_mask"][1:, 0].any()
