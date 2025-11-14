import math

import numpy as np

from motor_models import three_phase_waveform


def test_three_phase_waveform_balanced_set():
    t = np.linspace(0.0, 1 / 60, 360)
    va, vb, vc = three_phase_waveform(t, amplitude=1.2, frequency_hz=60.0, phase_offset_deg=10.0)
    assert math.isclose(np.max(va), 1.2, rel_tol=1e-3)
    assert math.isclose(np.max(vb), 1.2, rel_tol=1e-3)
    assert math.isclose(np.max(vc), 1.2, rel_tol=1e-3)
    space_vector_sum = va + vb + vc
    assert np.allclose(np.mean(space_vector_sum), 0.0, atol=1e-3)


def test_three_phase_waveform_includes_harmonic_component():
    t = np.linspace(0.0, 1 / 60, 360)
    base_va, *_ = three_phase_waveform(t, amplitude=1.0, frequency_hz=60.0, phase_offset_deg=0.0)
    harmonic_va, *_ = three_phase_waveform(
        t,
        amplitude=1.0,
        frequency_hz=60.0,
        phase_offset_deg=0.0,
        harmonic_order=5,
        harmonic_ratio=0.4,
    )
    assert not np.allclose(base_va, harmonic_va)
    assert math.isclose(np.std(harmonic_va - base_va), 0.4 * math.sqrt(0.5), rel_tol=0.2)
