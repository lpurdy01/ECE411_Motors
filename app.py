"""Dash application for visualizing ECE 411 synchronous machine concepts."""
from __future__ import annotations

import math
import os
import textwrap

import numpy as np
import plotly.graph_objects as go
from dash import Dash, Input, Output, dcc, html
from flask import abort, request

from motor_models import (
    FieldWeakeningParams,
    constant_power_torque,
    RoundRotorParams,
    SalientPoleParams,
    abc_to_alphabeta,
    abc_to_qd0,
    field_weakening_characteristics,
    field_weakening_flux,
    round_rotor_operating_point,
    salient_pole_torque_curves,
    three_phase_waveform,
)

APP_TITLE = "ECE 411 Motor Visualization"
TIME_VECTOR = np.linspace(0.0, 4.0 / 60.0, 600)
DELTA_RANGE = np.linspace(-math.pi, math.pi, 400)


def _phasor_trace(name: str, value: complex, color: str) -> go.Scatter:
    return go.Scatter(
        x=[0, value.real],
        y=[0, value.imag],
        mode="lines+markers",
        name=name,
        line=dict(color=color, width=4),
        marker=dict(color=color, size=8),
    )


def _phasor_layout(title: str) -> go.Layout:
    return go.Layout(
        title=title,
        xaxis=dict(title="α (real)", zeroline=True, range=[-1.6, 1.6]),
        yaxis=dict(title="β (imag)", zeroline=True, range=[-1.6, 1.6], scaleanchor="x"),
        showlegend=True,
        height=450,
    )


ACCESS_TOKEN = os.getenv("ACCESS_TOKEN")

app = Dash(__name__, suppress_callback_exceptions=True)
server = app.server

if ACCESS_TOKEN:

    @server.before_request
    def _check_token() -> None:
        token = request.args.get("token")
        if token != ACCESS_TOKEN:
            abort(403)

app.layout = html.Div(
    [
        html.H1(APP_TITLE),
        html.P(
            "Interactive playground for Clarke/Park transforms, synchronous "
            "machine torque, and field-weakening behavior."
        ),
        dcc.Tabs(
            id="main-tabs",
            value="tab-frames",
            children=[
                dcc.Tab(label="Frame Playground", value="tab-frames", children=[
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.H3("Signal controls"),
                                    dcc.Slider(
                                        id="frames-amplitude",
                                        min=0.2,
                                        max=1.8,
                                        value=1.0,
                                        step=0.05,
                                        marks=None,
                                        tooltip={"placement": "bottom"},
                                        updatemode="drag",
                                    ),
                                    html.Label("Peak phase voltage (pu)"),
                                    html.Br(),
                                    dcc.Slider(
                                        id="frames-frequency",
                                        min=10,
                                        max=120,
                                        value=60,
                                        step=1,
                                        marks={30: "30", 60: "60", 90: "90"},
                                    ),
                                    html.Label("Fundamental frequency (Hz)"),
                                    html.Br(),
                                    dcc.Slider(
                                        id="frames-phase",
                                        min=-180,
                                        max=180,
                                        value=0,
                                        step=5,
                                        marks={-180: "-180", 0: "0", 180: "180"},
                                    ),
                                    html.Label("Phase-a offset (deg)"),
                                    html.Br(),
                                    dcc.Dropdown(
                                        id="frames-harmonic-order",
                                        options=[
                                            {"label": "3rd", "value": 3},
                                            {"label": "5th", "value": 5},
                                            {"label": "7th", "value": 7},
                                        ],
                                        value=3,
                                        clearable=False,
                                        style={"marginBottom": "0.5rem"},
                                    ),
                                    dcc.Slider(
                                        id="frames-harmonic-ratio",
                                        min=0.0,
                                        max=0.6,
                                        value=0.0,
                                        step=0.05,
                                        marks={0: "0", 0.3: "0.3", 0.6: "0.6"},
                                    ),
                                    html.Label("Harmonic amplitude ratio"),
                                    html.Br(),
                                    dcc.Checklist(
                                        id="frames-sync",
                                        options=[{"label": "θₑ = ωₑt", "value": "sync"}],
                                        value=["sync"],
                                        style={"marginTop": "0.5rem"},
                                    ),
                                    dcc.Slider(
                                        id="frames-theta",
                                        min=0,
                                        max=360,
                                        value=0,
                                        step=5,
                                        marks={0: "0", 180: "180", 360: "360"},
                                    ),
                                    html.Label("Electrical angle θₑ (deg)"),
                                ],
                                className="tab-controls",
                            ),
                            html.Div(
                                [
                                    dcc.Graph(id="frames-abc-graph"),
                                    dcc.Graph(id="frames-alphabeta-graph"),
                                    dcc.Graph(id="frames-phasor-graph"),
                                    dcc.Graph(id="frames-qd-graph"),
                                    dcc.Markdown(
                                        id="frames-summary",
                                        mathjax=True,
                                        style={"backgroundColor": "#f8f9fa", "padding": "0.75rem"},
                                    ),
                                ],
                                className="tab-plots",
                            ),
                        ],
                        className="tab-layout",
                    )
                ]),
                dcc.Tab(label="Round-Rotor", value="tab-round", children=[
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.H3("Machine parameters"),
                                    html.Label("Rₛ (pu)"),
                                    dcc.Slider(id="round-rs", min=0.0, max=0.5, step=0.01, value=0.02),
                                    html.Label("Lₛ (pu)"),
                                    dcc.Slider(id="round-ls", min=0.1, max=2.0, step=0.05, value=0.6),
                                    html.Label("V₀ (pu)"),
                                    dcc.Slider(id="round-v0", min=0.3, max=1.5, step=0.05, value=1.0),
                                    html.Label("E₀ (pu)"),
                                    dcc.Slider(id="round-e0", min=0.3, max=1.5, step=0.05, value=0.9),
                                    html.Label("ωₑ (pu)"),
                                    dcc.Slider(id="round-we", min=0.3, max=2.0, step=0.05, value=1.0),
                                    html.Label("δ (deg)"),
                                    dcc.Slider(id="round-delta", min=-90, max=90, step=1, value=20),
                                ],
                                className="tab-controls",
                            ),
                            html.Div(
                                [
                                    dcc.Graph(id="round-phasor"),
                                    dcc.Markdown(
                                        id="round-summary",
                                        mathjax=True,
                                        style={"backgroundColor": "#f8f9fa", "padding": "0.75rem"},
                                    ),
                                ],
                                className="tab-plots",
                            ),
                        ],
                        className="tab-layout",
                    )
                ]),
                dcc.Tab(label="Field Weakening", value="tab-fw", children=[
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.H3("Limits"),
                                    html.Label("V₀,max (pu)"),
                                    dcc.Slider(id="fw-v0", min=0.5, max=1.5, step=0.05, value=1.0),
                                    html.Label("I₀,max (pu)"),
                                    dcc.Slider(id="fw-i0", min=0.2, max=2.0, step=0.05, value=1.2),
                                    html.Label("Lₛ (pu)"),
                                    dcc.Slider(id="fw-ls", min=0.1, max=2.0, step=0.05, value=0.6),
                                    html.Label("Rₛ (pu)"),
                                    dcc.Slider(id="fw-rs", min=0.0, max=0.3, step=0.01, value=0.0),
                                    html.Label("Operating ωₑ (pu)"),
                                    dcc.Slider(id="fw-speed", min=0.2, max=2.5, step=0.05, value=1.2),
                                ],
                                className="tab-controls",
                            ),
                            html.Div(
                                [
                                    dcc.Graph(id="fw-torque"),
                                    dcc.Graph(id="fw-power"),
                                    dcc.Graph(id="fw-flux"),
                                    dcc.Markdown(
                                        id="fw-summary",
                                        mathjax=True,
                                        style={"backgroundColor": "#f8f9fa", "padding": "0.75rem"},
                                    ),
                                ],
                                className="tab-plots",
                            ),
                        ],
                        className="tab-layout",
                    )
                ]),
                dcc.Tab(label="Salient Torque", value="tab-salient", children=[
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.H3("Salient parameters"),
                                    html.Label("L_d (pu)"),
                                    dcc.Slider(id="salient-ld", min=0.6, max=1.8, step=0.05, value=1.2),
                                    html.Label("L_q (pu)"),
                                    dcc.Slider(id="salient-lq", min=0.4, max=1.5, step=0.05, value=0.8),
                                    html.Label("E₀ (pu)"),
                                    dcc.Slider(id="salient-e0", min=0.3, max=1.5, step=0.05, value=1.0),
                                    html.Label("V₀ (pu)"),
                                    dcc.Slider(id="salient-v0", min=0.3, max=1.5, step=0.05, value=1.0),
                                    html.Label("ωₑ (pu)"),
                                    dcc.Slider(id="salient-we", min=0.3, max=2.0, step=0.05, value=1.0),
                                ],
                                className="tab-controls",
                            ),
                            html.Div(
                                [
                                    dcc.Graph(id="salient-graph"),
                                    dcc.Markdown(
                                        id="salient-summary",
                                        mathjax=True,
                                        style={"backgroundColor": "#f8f9fa", "padding": "0.75rem"},
                                    ),
                                ],
                                className="tab-plots",
                            ),
                        ],
                        className="tab-layout",
                    )
                ]),
            ],
        ),
    ],
    style={"padding": "1rem 2rem"},
)


@app.callback(
    Output("frames-abc-graph", "figure"),
    Output("frames-alphabeta-graph", "figure"),
    Output("frames-phasor-graph", "figure"),
    Output("frames-qd-graph", "figure"),
    Output("frames-summary", "children"),
    Input("frames-amplitude", "value"),
    Input("frames-frequency", "value"),
    Input("frames-phase", "value"),
    Input("frames-harmonic-order", "value"),
    Input("frames-harmonic-ratio", "value"),
    Input("frames-sync", "value"),
    Input("frames-theta", "value"),
)
def update_frames_tab(amplitude, frequency, phase, harmonic_order, harmonic_ratio, sync_value, theta_deg):
    va, vb, vc = three_phase_waveform(
        TIME_VECTOR,
        amplitude=amplitude,
        frequency_hz=frequency,
        phase_offset_deg=phase,
        harmonic_order=harmonic_order,
        harmonic_ratio=harmonic_ratio,
    )

    omega_e = 2 * math.pi * frequency
    if "sync" in (sync_value or []):
        theta = omega_e * TIME_VECTOR
    else:
        theta = np.full_like(TIME_VECTOR, math.radians(theta_deg))

    alpha, beta, v0 = abc_to_alphabeta(va, vb, vc)
    vq, vd = abc_to_qd0(va, vb, vc, theta)[:2]

    abc_fig = go.Figure()
    abc_fig.add_trace(go.Scatter(x=TIME_VECTOR * 1000, y=va, name="Va"))
    abc_fig.add_trace(go.Scatter(x=TIME_VECTOR * 1000, y=vb, name="Vb"))
    abc_fig.add_trace(go.Scatter(x=TIME_VECTOR * 1000, y=vc, name="Vc"))
    abc_fig.update_layout(title="Phase Voltages", xaxis_title="Time (ms)", yaxis_title="Voltage (pu)")

    ab_fig = go.Figure()
    ab_fig.add_trace(go.Scatter(x=TIME_VECTOR * 1000, y=alpha, name="α"))
    ab_fig.add_trace(go.Scatter(x=TIME_VECTOR * 1000, y=beta, name="β"))
    ab_fig.update_layout(title="αβ Components", xaxis_title="Time (ms)", yaxis_title="Magnitude (pu)")

    space_vector = alpha + 1j * beta
    phasor_fig = go.Figure()
    phasor_fig.add_trace(
        go.Scatter(
            x=np.real(space_vector),
            y=np.imag(space_vector),
            mode="lines",
            name="trajectory",
            line=dict(color="#7f7f7f"),
        )
    )
    phasor_fig.add_trace(_phasor_trace("instant", space_vector[-1], "#1f77b4"))
    phasor_fig.update_layout(_phasor_layout("Space Vector (αβ frame)"))

    qd_fig = go.Figure()
    qd_fig.add_trace(go.Scatter(x=TIME_VECTOR * 1000, y=vq, name="vq"))
    qd_fig.add_trace(go.Scatter(x=TIME_VECTOR * 1000, y=vd, name="vd"))
    qd_fig.update_layout(title="qd Components", xaxis_title="Time (ms)", yaxis_title="Magnitude (pu)")

    alpha_last = float(alpha[-1])
    beta_last = float(beta[-1])
    theta_last = float(np.atleast_1d(theta)[-1])
    vq_last = float(vq[-1])
    vd_last = float(vd[-1])
    space_rms = np.sqrt(np.mean(np.square(np.abs(space_vector))))

    summary = textwrap.dedent(
        fr"""
        **Clarke/Park snapshot**

        Balanced set? {'Yes' if np.allclose(v0, 0, atol=1e-3) else 'No'}  
        Space-vector rms ≈ {space_rms:.3f} pu

        $$\begin{{aligned}}
        v_q &= \alpha \cos\theta_e - \beta \sin\theta_e \\
            &= ({alpha_last:.3f})\cos({theta_last:.3f}) - ({beta_last:.3f})\sin({theta_last:.3f}) = {vq_last:.3f} \\
        v_d &= \alpha \sin\theta_e + \beta \cos\theta_e \\
            &= ({alpha_last:.3f})\sin({theta_last:.3f}) + ({beta_last:.3f})\cos({theta_last:.3f}) = {vd_last:.3f}
        \end{{aligned}}$$

        Mean values: $\overline{{v_q}} = {np.mean(vq):.3f}$, $\overline{{v_d}} = {np.mean(vd):.3f}$
        """
    ).strip()

    return abc_fig, ab_fig, phasor_fig, qd_fig, summary


@app.callback(
    Output("round-phasor", "figure"),
    Output("round-summary", "children"),
    Input("round-rs", "value"),
    Input("round-ls", "value"),
    Input("round-v0", "value"),
    Input("round-e0", "value"),
    Input("round-we", "value"),
    Input("round-delta", "value"),
)
def update_round_tab(rs, ls, v0, e0, w_e, delta_deg):
    params = RoundRotorParams(rs=rs, ls=ls, v0=v0, e0=e0, w_e=w_e)
    delta = math.radians(delta_deg)
    op = round_rotor_operating_point(params, delta)

    voltage = op["vq"] + 1j * op["vd"]
    emf = op["eq"] + 1j * op["ed"]
    current = op["iq"] + 1j * op["id"]
    drops_r = op["voltage_drop_res"]
    drops_l = op["voltage_drop_sync"]

    fig = go.Figure()
    fig.add_trace(_phasor_trace("Terminal V", voltage, "#1f77b4"))
    fig.add_trace(_phasor_trace("Internal E", emf, "#ff7f0e"))
    fig.add_trace(_phasor_trace("Stator I", current, "#2ca02c"))
    fig.add_trace(_phasor_trace("Rₛ I", drops_r, "#d62728"))
    fig.add_trace(_phasor_trace("jωLₛ I", drops_l, "#9467bd"))
    fig.update_layout(_phasor_layout("qd-Plane Phasors"))

    summary = textwrap.dedent(
        fr"""
        **Round-rotor steady state**

        $$\begin{{aligned}}
        V_q &= V_0 \cos\delta = {v0:.3f}\cos({delta:.3f}) = {op['vq']:.3f} \\
        V_d &= V_0 \sin\delta = {v0:.3f}\sin({delta:.3f}) = {op['vd']:.3f} \\
        P &= \tfrac{{3}}{{2}}(V_q I_q + V_d I_d) = \tfrac{{3}}{{2}}({op['vq']:.3f}\cdot{op['iq']:.3f} + {op['vd']:.3f}\cdot{op['id']:.3f}) = {op['power']:.3f} \\
        T_e &= \frac{{P}}{{\omega_e}} = \frac{{{op['power']:.3f}}}{{{w_e:.3f}}} = {op['torque']:.3f}
        \end{{aligned}}$$

        Power factor $\cos\varphi = {op['pf']:.3f}$ (ϕ = {math.degrees(op['pf_angle']):.1f}°)  
        $|I| = {op['i_mag']:.3f}$ pu
        """
    ).strip()

    return fig, summary


@app.callback(
    Output("fw-torque", "figure"),
    Output("fw-power", "figure"),
    Output("fw-flux", "figure"),
    Output("fw-summary", "children"),
    Input("fw-v0", "value"),
    Input("fw-i0", "value"),
    Input("fw-ls", "value"),
    Input("fw-rs", "value"),
    Input("fw-speed", "value"),
)
def update_field_weakening(v0_max, i0_max, ls, rs, selected_speed):
    params = FieldWeakeningParams(v0_max=v0_max, i0_max=i0_max, ls=ls, rs=rs)
    speed_grid = np.linspace(0.2, 2.5, 120)
    curves = field_weakening_characteristics(params, speed_grid)

    torque_fig = go.Figure()
    torque_fig.add_trace(
        go.Scatter(x=curves["speed"], y=curves["torque"], mode="lines", name="Tₑ,max")
    )
    torque_fig.add_vline(x=selected_speed, line_dash="dash", line_color="#d62728")
    torque_fig.update_layout(title="Torque Limit vs Speed", xaxis_title="ωₑ (pu)", yaxis_title="Torque (pu)")

    power_fig = go.Figure()
    power_fig.add_trace(
        go.Scatter(x=curves["speed"], y=curves["power"], mode="lines", name="Pₑ,max", line_color="#ff7f0e")
    )
    power_fig.add_vline(x=selected_speed, line_dash="dash", line_color="#d62728")
    power_fig.update_layout(title="Power Limit vs Speed", xaxis_title="ωₑ (pu)", yaxis_title="Power (pu)")

    flux_fig = go.Figure()
    flux_fig.add_trace(
        go.Scatter(x=curves["speed"], y=curves["flux"], mode="lines", name="|Λ₀|")
    )
    flux_fig.add_vline(x=selected_speed, line_dash="dash", line_color="#d62728")
    flux_fig.update_layout(title="Flux Linkage vs Speed", xaxis_title="ωₑ (pu)", yaxis_title="|Λ₀| (pu)")

    flux_selected = field_weakening_flux(v0_max, i0_max, ls, selected_speed)
    torque_selected, power_selected = constant_power_torque(v0_max, i0_max, selected_speed) if rs == 0 else (
        curves["torque"][np.argmin(np.abs(curves["speed"] - selected_speed))],
        curves["power"][np.argmin(np.abs(curves["speed"] - selected_speed))],
    )

    if rs == 0:
        summary = textwrap.dedent(
            fr"""
            **Field-weakening snapshot**

            $$\Lambda_0 = \sqrt{{\left(\frac{{V_0}}{{\omega_e}}\right)^2 + (L_s I_0)^2}} = \sqrt{{\left(\frac{{{v0_max:.3f}}}{{{selected_speed:.3f}}}\right)^2 + ({ls:.3f}\cdot{i0_max:.3f})^2}} = {flux_selected:.3f}$$

            $$T_{{e,\max}} = \frac{{V_0 I_0}}{{\omega_e}} = \frac{{{v0_max:.3f}\cdot{i0_max:.3f}}}{{{selected_speed:.3f}}} = {torque_selected:.3f}$$

            $$P_{{\max}} = V_0 I_0 = {v0_max:.3f}\cdot{i0_max:.3f} = {power_selected:.3f}$$
            """
        ).strip()
    else:
        zs = math.hypot(rs, selected_speed * ls)
        zs_sq = zs**2
        e0 = selected_speed * flux_selected
        summary = textwrap.dedent(
            fr"""
            **Field-weakening snapshot (Rₛ ≠ 0)**

            $$Z_s = \sqrt{{R_s^2 + (\omega_e L_s)^2}} = \sqrt{{{rs:.3f}^2 + ({selected_speed:.3f}\cdot{ls:.3f})^2}} = {zs:.3f}$$

            $$T_{{e,\max}} = \frac{{V_0 E_0}}{{\omega_e Z_s}} - \frac{{E_0^2 R_s}}{{\omega_e Z_s^2}} = \frac{{{v0_max:.3f}\cdot{e0:.3f}}}{{{selected_speed:.3f}\cdot{zs:.3f}}} - \frac{{{e0:.3f}^2\cdot{rs:.3f}}}{{{selected_speed:.3f}\cdot{zs_sq:.3f}}} = {torque_selected:.3f}$$

            $$P_{{\max}} = T_{{e,\max}}\, \omega_e = {torque_selected:.3f}\cdot{selected_speed:.3f} = {power_selected:.3f}$$
            """
        ).strip()

    return torque_fig, power_fig, flux_fig, summary


@app.callback(
    Output("salient-graph", "figure"),
    Output("salient-summary", "children"),
    Input("salient-ld", "value"),
    Input("salient-lq", "value"),
    Input("salient-e0", "value"),
    Input("salient-v0", "value"),
    Input("salient-we", "value"),
)
def update_salient_tab(ld, lq, e0, v0, w_e):
    params = SalientPoleParams(ld=ld, lq=lq, e0=e0, v0=v0, w_e=w_e)
    t_field, t_rel, t_total = salient_pole_torque_curves(params, DELTA_RANGE)

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=np.degrees(DELTA_RANGE), y=t_field, name="Field torque"))
    fig.add_trace(go.Scatter(x=np.degrees(DELTA_RANGE), y=t_rel, name="Reluctance torque"))
    fig.add_trace(go.Scatter(x=np.degrees(DELTA_RANGE), y=t_total, name="Total torque", line=dict(width=3)))
    fig.update_layout(
        title="Torque vs Load Angle",
        xaxis_title="δ (deg)",
        yaxis_title="Torque (pu)",
        legend=dict(orientation="h"),
    )

    field_max = np.max(np.abs(t_field))
    rel_max = np.max(np.abs(t_rel))
    total_max = np.max(np.abs(t_total))
    saliency = ld / lq if lq != 0 else float("inf")

    summary = textwrap.dedent(
        fr"""
        **Salient-pole torque**

        $$T_{{field}}(\delta) = -\frac{{E_0 V_0}}{{\omega_e L_d}} \sin\delta = -\frac{{{e0:.3f}\cdot{v0:.3f}}}{{{w_e:.3f}\cdot{ld:.3f}}} \sin\delta$$

        $$T_{{rel}}(\delta) = -\frac{{V_0^2 (L_d - L_q)}}{{2\, \omega_e L_d L_q}} \sin(2\delta) = -\frac{{{v0:.3f}^2 ({ld:.3f}-{lq:.3f})}}{{2\cdot{w_e:.3f}\cdot{ld:.3f}\cdot{lq:.3f}}} \sin(2\delta)$$

        Peak magnitudes: $|T_{{field}}|_{{\max}} \approx {field_max:.3f}$ pu, $|T_{{rel}}|_{{\max}} \approx {rel_max:.3f}$ pu, $|T_{{total}}|_{{\max}} \approx {total_max:.3f}$ pu.  
        Saliency ratio $\frac{{L_d}}{{L_q}} = {saliency:.3f}$
        """
    ).strip()

    return fig, summary


if __name__ == "__main__":
    app.run(debug=True)
