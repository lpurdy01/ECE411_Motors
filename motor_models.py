"""Core mathematical models and helper utilities for the ECE 411 motor
visualization dash application.

This module intentionally contains no Dash imports so that it can be unit
tested in isolation.  All functions operate in per-unit wherever practical and
follow the Lipo–Novotny convention for Clarke and Park transforms.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable, Tuple

import numpy as np


SQRT3 = math.sqrt(3.0)
TWO_THIRDS = 2.0 / 3.0
THREE_HALVES = 1.5


@dataclass(frozen=True)
class MachineParams:
    """Physical parameters for the synchronous machine and inverter.

    Defaults mirror ``machine_params.m`` from the MATLAB studies so that the
    Dash UI and tests start from the same baseline.  All values are expressed
    in SI units:

    * ``Lm`` magnetizing inductance [H]
    * ``Ls`` stator synchronous inductance (q = d) [H]
    * ``Rs`` stator resistance [Ω]
    * ``Rf`` field resistance [Ω]
    * ``p`` pole count
    * ``wb`` electrical base speed [rad/s]
    * ``wmb`` mechanical base speed [rad/s]
    * ``Vdc`` dc bus voltage [V]
    * ``I0_max`` peak stator current limit [A]
    * ``If_max`` dc field current limit [A]
    """

    Lm: float = 17.5e-3
    Ls: float = 0.6e-3
    Rs: float = 0.0155
    Rf: float = 26.5
    p: int = 8
    wb: float = 2 * math.pi * 267
    wmb: float = 419.0
    Vdc: float = 550.0
    I0_max: float = 300.0
    If_max: float = 7.0


@dataclass(frozen=True)
class PUBases:
    """Per-unit bases following the Part 1 MATLAB setup.

    The Park-frame voltage base is line-neutral peak so that the modulation
    index relationship uses ``M0 = 2/√3 · V0_pu``.
    """

    VB: float  # LN peak base voltage [V]
    PB: float  # three-phase base power [W]
    wB: float  # electrical base speed [rad/s]
    IB: float  # phase current peak base [A]
    ZB: float  # impedance base [Ω]
    LB: float  # inductance base [H]
    wmB: float  # mechanical base speed [rad/s]
    TB: float  # torque base [N·m]


@dataclass(frozen=True)
class RoundRotorParams:
    """Parameters describing a round-rotor synchronous machine in per-unit."""

    rs: float
    ls: float
    v0: float
    e0: float
    w_e: float


@dataclass(frozen=True)
class SalientPoleParams:
    """Parameters describing a salient-pole synchronous machine in per-unit."""

    ld: float
    lq: float
    e0: float
    v0: float
    w_e: float


@dataclass(frozen=True)
class FieldWeakeningParams:
    """Parameters for field weakening studies."""

    v0_max: float
    i0_max: float
    ls: float
    rs: float = 0.0


# ---------------------------------------------------------------------------
# Clarke and Park transforms
# ---------------------------------------------------------------------------

def abc_to_alphabeta(va: np.ndarray, vb: np.ndarray, vc: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return alpha, beta, v0 using the Lipo–Novotny Clarke transform."""

    alpha = TWO_THIRDS * (va - 0.5 * vb - 0.5 * vc)
    beta = TWO_THIRDS * (SQRT3 / 2.0) * (vc - vb)
    v0 = (va + vb + vc) / 3.0
    return alpha, beta, v0


def alphabeta_to_qd(alpha: np.ndarray, beta: np.ndarray, theta_e: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Convert alpha/beta components to q/d components."""

    c = np.cos(theta_e)
    s = np.sin(theta_e)
    vq = alpha * c - beta * s
    vd = alpha * s + beta * c
    return vq, vd


def qd_to_alphabeta(vq: np.ndarray, vd: np.ndarray, theta_e: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Inverse Park transform."""

    c = np.cos(theta_e)
    s = np.sin(theta_e)
    alpha = vq * c + vd * s
    beta = -vq * s + vd * c
    return alpha, beta


def alphabeta_to_abc(alpha: np.ndarray, beta: np.ndarray, v0: np.ndarray | float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Inverse Clarke transform returning phase voltages."""

    va = alpha + v0
    vb = -0.5 * alpha - (SQRT3 / 2.0) * beta + v0
    vc = -0.5 * alpha + (SQRT3 / 2.0) * beta + v0
    return va, vb, vc


def abc_to_qd0(va: np.ndarray, vb: np.ndarray, vc: np.ndarray, theta_e: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convenience wrapper performing Clarke then Park transform."""

    alpha, beta, v0 = abc_to_alphabeta(va, vb, vc)
    vq, vd = alphabeta_to_qd(alpha, beta, theta_e)
    return vq, vd, v0


# ---------------------------------------------------------------------------
# Waveform helpers
# ---------------------------------------------------------------------------

def three_phase_waveform(
    t: np.ndarray,
    amplitude: float = 1.0,
    frequency_hz: float = 60.0,
    phase_offset_deg: float = 0.0,
    harmonic_order: int = 1,
    harmonic_ratio: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Generate three-phase waveforms with optional harmonic content."""

    omega = 2.0 * math.pi * frequency_hz
    phase_offset = math.radians(phase_offset_deg)

    base = amplitude * np.cos(omega * t + phase_offset)
    phase_shift = 2.0 * math.pi / 3.0
    va = base
    vb = amplitude * np.cos(omega * t + phase_offset - phase_shift)
    vc = amplitude * np.cos(omega * t + phase_offset + phase_shift)

    if harmonic_ratio != 0.0:
        h = float(harmonic_order)
        harmonic = harmonic_ratio * amplitude
        va += harmonic * np.cos(h * omega * t + phase_offset)
        vb += harmonic * np.cos(h * omega * t + phase_offset - phase_shift * h)
        vc += harmonic * np.cos(h * omega * t + phase_offset + phase_shift * h)

    return va, vb, vc


def default_machine_params() -> MachineParams:
    """Return the default :class:`MachineParams` used across the MATLAB demos."""

    return MachineParams()


def pu_bases_from_machine(params: MachineParams) -> PUBases:
    """Construct per-unit bases consistent with the MATLAB notebooks.

    Formulas (LN denotes line-to-neutral):

    * Voltage base (LN peak): ``V_B = V_dc/√3``
    * Power base: ``P_B = 1.5 · V_B · I_0,max``
    * Current base (LN peak): ``I_B = 2 P_B / (3 V_B)``
    * Impedance base: ``Z_B = V_B / I_B``
    * Inductance base: ``L_B = Z_B / ω_B``
    * Mechanical base speed: ``ω_m,B = 2 ω_B / p``
    * Torque base: ``T_B = (p/2)(P_B/ω_B)``
    """

    vb_ln_peak = params.Vdc / math.sqrt(3.0)
    vb_ln_rms = vb_ln_peak / math.sqrt(2.0)
    pb_3ph = 1.5 * vb_ln_peak * params.I0_max

    vb = math.sqrt(2.0) * vb_ln_rms
    pb = pb_3ph
    wb = params.wb
    ib = (2.0 * pb) / (3.0 * vb)
    zb = vb / ib
    lb = zb / wb
    wmb = 2.0 * wb / params.p
    tb = (params.p / 2.0) * (pb / wb)

    return PUBases(VB=vb, PB=pb, wB=wb, IB=ib, ZB=zb, LB=lb, wmB=wmb, TB=tb)


# ---------------------------------------------------------------------------
# Round-rotor synchronous machine
# ---------------------------------------------------------------------------

def round_rotor_operating_point(params: RoundRotorParams, delta: float) -> dict:
    """Compute the operating point quantities for a round-rotor machine."""

    vq = params.v0 * math.cos(delta)
    vd = params.v0 * math.sin(delta)
    eq = params.e0
    ed = 0.0

    wls = params.w_e * params.ls
    matrix = np.array([[params.rs, wls], [-wls, params.rs]], dtype=float)
    rhs = np.array([vq - eq, vd - ed], dtype=float)
    iq, id_ = np.linalg.solve(matrix, rhs)

    i_mag = math.hypot(iq, id_)
    power = THREE_HALVES * (vq * iq + vd * id_)
    torque = power / params.w_e if params.w_e != 0 else float("nan")

    voltage_drop_res = params.rs * (iq + 1j * id_)
    voltage_drop_sync = 1j * params.w_e * params.ls * (iq + 1j * id_)

    pf_angle = math.atan2(id_, iq)
    pf = math.cos(pf_angle)

    return {
        "vq": vq,
        "vd": vd,
        "iq": iq,
        "id": id_,
        "i_mag": i_mag,
        "power": power,
        "torque": torque,
        "pf_angle": pf_angle,
        "pf": pf,
        "voltage_drop_res": voltage_drop_res,
        "voltage_drop_sync": voltage_drop_sync,
        "eq": eq,
        "ed": ed,
}


def sm_required_voltage_park(
    mechanical_omega_m: np.ndarray | float,
    field_current: float,
    iq: np.ndarray | float,
    id_: np.ndarray | float,
    params: MachineParams,
) -> dict:
    """Stator voltage required for a given operating point (SI units).

    Park-frame KVL: ``V_qd = E_qd + (R_s + j ω_e L_s) I_qd`` with Lipo–Novotny
    axes (``q = Re{·}``, ``d = -Im{·}``).  The helper returns a dictionary with
    ``V``, components ``Vq``/``Vd``, magnitude ``V0``, electrical speed
    ``w_e``, stator impedance ``Zs``, and internal EMF ``E``.
    """

    w_e = (params.p / 2.0) * np.asarray(mechanical_omega_m, dtype=float)
    zs = params.Rs + 1j * w_e * params.Ls
    e_internal = (w_e * params.Lm * field_current) + 1j * 0.0

    current = np.asarray(iq, dtype=float) + 1j * (-np.asarray(id_, dtype=float))
    v_phasor = e_internal + zs * current

    return {
        "V": v_phasor,
        "Vq": np.real(v_phasor),
        "Vd": -np.imag(v_phasor),
        "V0": np.abs(v_phasor),
        "w_e": w_e,
        "Zs": zs,
        "E": e_internal,
    }


def sm_required_voltage_park_pu(
    mechanical_omega_m: np.ndarray | float,
    field_current: float,
    iq: np.ndarray | float,
    id_: np.ndarray | float,
    params: MachineParams,
    bases: PUBases,
) -> dict:
    """Per-unit Park-frame voltage calculation.

    The computation mirrors :func:`sm_required_voltage_park` and rescales each
    quantity by the provided per-unit bases.  ``V0_pu`` uses the line-to-neutral
    peak base consistent with the MATLAB scripts.
    """

    w_e = (params.p / 2.0) * np.asarray(mechanical_omega_m, dtype=float)
    w_pu = w_e / bases.wB

    rs_pu = params.Rs / bases.ZB
    xs_base_pu = (bases.wB * params.Ls) / bases.ZB

    zs_pu = rs_pu + 1j * (w_pu * xs_base_pu)
    e_pu = (w_e * params.Lm * field_current) / bases.VB

    current_pu = (np.asarray(iq, dtype=float) + 1j * (-np.asarray(id_, dtype=float))) / bases.IB
    v_pu = e_pu + zs_pu * current_pu

    return {
        "V_pu": v_pu,
        "Vq_pu": np.real(v_pu),
        "Vd_pu": -np.imag(v_pu),
        "V0_pu": np.abs(v_pu),
        "w_e": w_e,
        "Zs_pu": zs_pu,
        "E_pu": e_pu,
    }


# ---------------------------------------------------------------------------
# Salient-pole torque decomposition
# ---------------------------------------------------------------------------

def salient_pole_torque_curves(params: SalientPoleParams, delta: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return field, reluctance, and total torque for an array of load angles."""

    t_field = -(params.e0 * params.v0) / (params.w_e * params.ld) * np.sin(delta)
    t_rel = -(params.v0**2 * (params.ld - params.lq)) / (2.0 * params.w_e * params.ld * params.lq) * np.sin(2.0 * delta)
    t_total = t_field + t_rel
    return t_field, t_rel, t_total


# ---------------------------------------------------------------------------
# Field-weakening helpers
# ---------------------------------------------------------------------------

def field_weakening_flux(v0: float, i0: float, ls: float, w_e: float) -> float:
    """Compute the required flux linkage magnitude under field weakening."""

    return math.sqrt((v0 / w_e) ** 2 + (ls * i0) ** 2)


def constant_power_torque(v0: float, i0: float, w_e: float) -> Tuple[float, float]:
    """Return maximum torque and power in the constant-power region."""

    p_max = v0 * i0
    te_max = p_max / w_e if w_e != 0 else float("inf")
    return te_max, p_max


def field_weakening_characteristics(
    params: FieldWeakeningParams,
    speed_grid: Iterable[float],
) -> dict:
    """Compute torque, power, and flux characteristics across speed."""

    speeds = np.asarray(list(speed_grid), dtype=float)
    torque = np.zeros_like(speeds)
    power = np.zeros_like(speeds)
    flux = np.zeros_like(speeds)

    for idx, w_e in enumerate(speeds):
        if w_e <= 0:
            flux[idx] = float("nan")
            torque[idx] = 0.0
            power[idx] = 0.0
            continue

        flux_val = field_weakening_flux(params.v0_max, params.i0_max, params.ls, w_e)
        flux[idx] = flux_val

        if params.rs == 0:
            torque[idx], power[idx] = constant_power_torque(params.v0_max, params.i0_max, w_e)
        else:
            zs = math.hypot(params.rs, w_e * params.ls)
            e0 = w_e * flux_val
            torque[idx] = (params.v0_max * e0) / (w_e * zs) - (e0**2 * params.rs) / (w_e * (zs**2))
            power[idx] = torque[idx] * w_e

    return {"speed": speeds, "torque": torque, "power": power, "flux": flux}


def inverter_voltage_vs_speed(
    params: MachineParams,
    bases: PUBases,
    field_currents: Iterable[float],
    omega_m: np.ndarray,
) -> dict:
    """Voltage requirement sweep vs mechanical speed for multiple ``I_f``.

    Replicates MATLAB Part 1 by holding ``Iq = I0_max`` and ``Id = 0`` while
    evaluating :math:`V_0 = |E_q + (R_s + j\omega_e L_s) I_qd|`.
    """

    v0_lim = params.Vdc / math.sqrt(3.0)
    v0_lim_pu = v0_lim / bases.VB

    traces = []
    for ifield in field_currents:
        result = sm_required_voltage_park_pu(omega_m, ifield, params.I0_max, 0.0, params, bases)
        v0_pu = np.asarray(result["V0_pu"], dtype=float)
        cross_idx = int(np.argmax(v0_pu >= v0_lim_pu)) if np.any(v0_pu >= v0_lim_pu) else None
        traces.append(
            {
                "If": ifield,
                "omega_m": omega_m,
                "V0": v0_pu * bases.VB,
                "V0_pu": v0_pu,
                "cross_idx": cross_idx,
                "cross_speed": omega_m[cross_idx] if cross_idx is not None else None,
            }
        )

    return {"limit": v0_lim, "limit_pu": v0_lim_pu, "traces": traces}


def mechanical_power_map(
    field_current: float,
    torque_values: np.ndarray,
    omega_m_values: np.ndarray,
    params: MachineParams,
    bases: PUBases,
) -> dict:
    """Compute mechanical power contours (Part 2) for a fixed field current.

    The torque requirement uses ``I_q = (T_e · ω_e) / E_q`` with ``I_d = 0`` and
    ``E_q = ω_e L_m I_f``.  Feasible regions are bounded by current and voltage
    limits.
    """

    omega_m_grid, torque_grid = np.meshgrid(omega_m_values, torque_values)
    omega_e_grid = (params.p / 2.0) * omega_m_grid
    omega_e_pu = omega_e_grid / bases.wB

    eq_pu = (omega_e_grid * params.Lm * field_current) / bases.VB
    eq_pu[np.abs(eq_pu) < 1e-9] = np.nan

    iq_pu = (torque_grid / bases.TB * omega_e_pu) / eq_pu
    iq_pu[~np.isfinite(iq_pu)] = np.nan

    v_pu = sm_required_voltage_park_pu(omega_m_grid, field_current, iq_pu * bases.IB, 0.0, params, bases)
    v0_pu = np.reshape(v_pu["V0_pu"], omega_m_grid.shape)

    i_mag_pu = np.abs(iq_pu)
    current_mask = i_mag_pu <= (params.I0_max / bases.IB)
    voltage_mask = v0_pu <= (params.Vdc / math.sqrt(3.0) / bases.VB)
    feasible_mask = current_mask & voltage_mask & np.isfinite(iq_pu)

    p_mech_pu = np.real(eq_pu) * iq_pu
    p_mech_w = p_mech_pu * bases.PB
    p_mech_w[~feasible_mask] = np.nan

    return {
        "omega_m": omega_m_grid,
        "torque": torque_grid,
        "power_w": p_mech_w,
        "voltage_mask": voltage_mask,
        "current_mask": current_mask,
        "feasible_mask": feasible_mask,
    }


def modulation_index_map(
    field_current: float,
    torque_values: np.ndarray,
    omega_m_values: np.ndarray,
    params: MachineParams,
    bases: PUBases,
) -> dict:
    """Compute modulation index contours ``M0 = 2/√3 · V0_pu`` for Part 3.

    Values outside the simultaneous voltage/current limits are set to NaN so
    contour plots clip to the feasible region. Returned mask and headroom
    arrays allow the UI to render constraint overlays without recomputing.
    """

    omega_m_grid, torque_grid = np.meshgrid(omega_m_values, torque_values)
    omega_e_grid = (params.p / 2.0) * omega_m_grid
    eq_pu = (omega_e_grid * params.Lm * field_current) / bases.VB
    eq_pu[np.abs(eq_pu) < 1e-9] = np.nan

    omega_e_pu = omega_e_grid / bases.wB
    torque_pu = torque_grid / bases.TB
    iq_pu = (torque_pu * omega_e_pu) / eq_pu
    iq_pu[~np.isfinite(iq_pu)] = np.nan

    v_pu = sm_required_voltage_park_pu(omega_m_grid, field_current, iq_pu * bases.IB, 0.0, params, bases)
    v0_pu = np.reshape(v_pu["V0_pu"], omega_m_grid.shape)
    m0 = (2.0 / math.sqrt(3.0)) * v0_pu
    i_mag_pu = np.abs(iq_pu)
    current_ratio = i_mag_pu / (params.I0_max / bases.IB)
    voltage_ratio = v0_pu / (params.Vdc / math.sqrt(3.0) / bases.VB)

    current_mask = i_mag_pu <= (params.I0_max / bases.IB)
    voltage_mask = v0_pu <= (params.Vdc / math.sqrt(3.0) / bases.VB)
    feasible_mask = current_mask & voltage_mask & np.isfinite(iq_pu)

    m0_masked = np.where(feasible_mask, m0, np.nan)

    return {
        "omega_m": omega_m_grid,
        "torque": torque_grid,
        "m0": m0,
        "m0_masked": m0_masked,
        "current_mask": current_mask,
        "voltage_mask": voltage_mask,
        "feasible_mask": feasible_mask,
        "current_ratio": current_ratio,
        "voltage_ratio": voltage_ratio,
    }


def loss_curves(params: MachineParams) -> dict:
    """Return 1D loss curves for field, core, and stator/inverter blocks."""

    if_axis = np.linspace(0.0, params.If_max, 400)
    pf_axis = (if_axis**2) * params.Rf

    e0_rated = params.wb * params.Lm * params.If_max
    e0_axis = np.linspace(0.0, e0_rated, 400)
    pcore_axis = 500.0 * (e0_axis / e0_rated)

    i0_axis = np.linspace(0.0, params.I0_max, 400)
    ps_axis = 0.5 * (i0_axis**2) * params.Rs
    pinv_axis = 0.015 * params.Vdc * i0_axis
    ps_total_axis = ps_axis + pinv_axis

    return {
        "If": if_axis,
        "Pf": pf_axis,
        "E0": e0_axis,
        "Pcore": pcore_axis,
        "I0": i0_axis,
        "Ps": ps_axis,
        "Pinv": pinv_axis,
        "Ps_total": ps_total_axis,
        "E0_rated": e0_rated,
    }


def loss_surface(
    field_current: float,
    torque_values: np.ndarray,
    omega_m_values: np.ndarray,
    params: MachineParams,
) -> dict:
    """Compute total loss surface for a fixed field current (extra credit plot)."""

    omega_m_grid, torque_grid = np.meshgrid(omega_m_values, torque_values)
    flux_linkage = params.Lm * field_current

    omega_e_grid = (params.p / 2.0) * omega_m_grid
    eq_grid = omega_e_grid * flux_linkage

    iq_grid = (4.0 / (3.0 * params.p)) * (torque_grid / flux_linkage)
    i_mag_grid = np.abs(iq_grid)

    vq_grid = iq_grid * params.Rs + eq_grid
    vd_grid = -omega_e_grid * params.Ls * iq_grid
    v0_grid = np.sqrt(vq_grid**2 + vd_grid**2)

    current_mask = i_mag_grid <= params.I0_max
    voltage_mask = v0_grid <= (params.Vdc / math.sqrt(3.0))
    feasible_mask = current_mask & voltage_mask

    pf_scalar = (field_current**2) * params.Rf
    e0_grid = np.abs(eq_grid)
    pcore_grid = 500.0 * (e0_grid / (params.wb * params.Lm * params.If_max))
    ps_grid = 0.5 * (i_mag_grid**2) * params.Rs
    pinv_grid = 0.015 * params.Vdc * i_mag_grid

    p_total_grid = pf_scalar + pcore_grid + ps_grid + pinv_grid
    p_total_grid[~feasible_mask] = np.nan

    return {
        "omega_m": omega_m_grid,
        "torque": torque_grid,
        "loss_kw": p_total_grid / 1e3,
        "feasible_mask": feasible_mask,
    }


__all__ = [
    "MachineParams",
    "PUBases",
    "RoundRotorParams",
    "SalientPoleParams",
    "FieldWeakeningParams",
    "abc_to_alphabeta",
    "alphabeta_to_qd",
    "qd_to_alphabeta",
    "alphabeta_to_abc",
    "abc_to_qd0",
    "three_phase_waveform",
    "default_machine_params",
    "pu_bases_from_machine",
    "round_rotor_operating_point",
    "salient_pole_torque_curves",
    "field_weakening_flux",
    "constant_power_torque",
    "field_weakening_characteristics",
    "sm_required_voltage_park",
    "sm_required_voltage_park_pu",
    "inverter_voltage_vs_speed",
    "mechanical_power_map",
    "modulation_index_map",
    "loss_curves",
    "loss_surface",
]
