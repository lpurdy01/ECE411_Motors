import math

from app import (
    update_field_weakening,
    update_frames_tab,
    update_round_tab,
    update_salient_tab,
)


def test_update_frames_tab_returns_consistent_figures():
    outputs = update_frames_tab(1.0, 60, 0, 3, 0.0, ["sync"], 0)
    abc_fig, ab_fig, phasor_fig, qd_fig, summary = outputs
    assert all(fig.data for fig in (abc_fig, ab_fig, qd_fig))
    assert "Balanced set? Yes" in summary
    assert "v_q" in summary


def test_update_frames_tab_manual_theta_changes_summary():
    _, _, _, _, summary_sync = update_frames_tab(1.0, 60, 0, 3, 0.0, ["sync"], 0)
    _, _, _, _, summary_manual = update_frames_tab(1.0, 60, 0, 3, 0.0, [], 90)
    assert summary_sync != summary_manual


def test_update_round_tab_summary_contains_expected_values():
    fig, summary = update_round_tab(0.1, 0.5, 1.0, 0.9, 1.0, 15)
    assert fig.data
    assert "Round-rotor steady state" in summary
    assert "V_q" in summary and "V_d" in summary


def test_update_field_weakening_res_zero_matches_closed_form():
    torque_fig, power_fig, flux_fig, summary = update_field_weakening(1.0, 0.8, 0.5, 0.0, 1.2)
    assert torque_fig.data and power_fig.data and flux_fig.data
    expected_flux = math.sqrt((1.0 / 1.2) ** 2 + (0.5 * 0.8) ** 2)
    assert f"{expected_flux:.3f}" in summary


def test_update_salient_tab_summary_mentions_peaks():
    fig, summary = update_salient_tab(1.1, 0.8, 1.0, 0.95, 1.0)
    assert fig.data
    assert "Salient-pole torque" in summary
    assert "Saliency ratio" in summary
    assert any(trace.name == "Total torque" for trace in fig.data)
