import math

import math
import numpy as np

from motor_models import (
    default_machine_params,
    efficiency_map,
    inverter_voltage_vs_speed,
    mechanical_power_map,
    modulation_index_map,
    pu_bases_from_machine,
    sm_required_voltage_park,
    sm_required_voltage_park_pu,
)
import plotly.graph_objects as go

from motor_tab import (
    THREE_D_SENTINEL,
    _headroom_from_mod_map,
    compute_motor_maps,
    update_motor_plots,
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


def test_modulation_map_clips_to_feasible_region():
    params = default_machine_params()
    bases = pu_bases_from_machine(params)
    torque_axis = np.linspace(0.0, 250.0, 10)
    omega_axis = np.linspace(0.0, 200.0, 10)

    mod_map = modulation_index_map(params.If_max, torque_axis, omega_axis, params, bases)
    # Regions violating voltage should be masked in m0_masked
    assert np.isnan(mod_map["m0_masked"][0, -1]) or mod_map["voltage_mask"][0, -1]

    # Pick a small but nonzero operating point to avoid the ω_e = 0 singularity
    feasible_point = mod_map["m0_masked"][1, 1]
    assert np.isfinite(feasible_point)


def test_efficiency_map_is_bounded_and_masked():
    params = default_machine_params()
    bases = pu_bases_from_machine(params)
    torque_axis = np.linspace(0.0, 200.0, 15)
    omega_axis = np.linspace(0.0, 300.0, 15)

    eff_map = efficiency_map(params.If_max, torque_axis, omega_axis, params, bases)
    assert np.nanmax(eff_map["efficiency"]) <= 1.0 + 1e-6
    assert np.isnan(eff_map["efficiency"][0, 0])


def test_headroom_helper_caps_values_and_masks():
    mod_map = {
        "current_ratio": np.array([[0.1, 0.0], [1.0, 2.0]]),
        "voltage_ratio": np.array([[0.1, 0.0], [0.5, 0.2]]),
        "feasible_mask": np.array([[True, False], [True, True]]),
    }

    headroom = _headroom_from_mod_map(mod_map)
    assert np.nanmax(headroom) <= 5.0
    assert np.isnan(headroom[0, 1])


def test_compute_motor_maps_includes_surface_option():
    params = default_machine_params()
    payload, options, selected, *_ = compute_motor_maps(
        0,
        params.Vdc,
        params.Ls,
        params.Lm,
        params.Rs,
        params.Rf,
        params.p,
        params.wb,
        params.I0_max,
        params.If_max,
        "",
        ["auto"],
        params.wmb,
        200,
        80,
        80,
    )

    assert options[0]["value"] == THREE_D_SENTINEL
    assert selected != THREE_D_SENTINEL
    assert payload["per_if"]


def test_surface_lookup_handles_integer_field_case_keys():
    params = default_machine_params()
    payload, *_ = compute_motor_maps(
        0,
        params.Vdc,
        params.Ls,
        params.Lm,
        params.Rs,
        params.Rf,
        params.p,
        params.wb,
        params.I0_max,
        params.If_max,
        "7,4,2",
        ["auto"],
        params.wmb,
        200,
        60,
        60,
    )

    # Simulate JSON round-trip that strips trailing decimals ("7.0" -> "7")
    payload["per_if"] = {key.rstrip("0").rstrip("."): val for key, val in payload["per_if"].items()}

    figures = update_motor_plots(THREE_D_SENTINEL, payload)
    surface_figs = figures[10:]

    assert all(isinstance(fig, go.Figure) for fig in surface_figs)
