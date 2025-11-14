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


__all__ = [
    "RoundRotorParams",
    "SalientPoleParams",
    "FieldWeakeningParams",
    "abc_to_alphabeta",
    "alphabeta_to_qd",
    "qd_to_alphabeta",
    "alphabeta_to_abc",
    "abc_to_qd0",
    "three_phase_waveform",
    "round_rotor_operating_point",
    "salient_pole_torque_curves",
    "field_weakening_flux",
    "constant_power_torque",
    "field_weakening_characteristics",
]
