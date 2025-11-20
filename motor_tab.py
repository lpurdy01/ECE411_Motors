"""Motor calculator layout and callbacks for the Dash app."""
from __future__ import annotations

import io
import json
import logging
import math
import time
import zipfile
from dataclasses import asdict

import numpy as np
import plotly.graph_objects as go
from dash import Input, Output, State, callback, dcc, html, no_update

from motor_models import (
    MachineParams,
    PUBases,
    default_machine_params,
    efficiency_map,
    inverter_voltage_vs_speed,
    loss_curves,
    loss_surface,
    mechanical_power_map,
    modulation_index_map,
    pu_bases_from_machine,
)

DEFAULT_MACHINE = default_machine_params()
DEFAULT_BASES = pu_bases_from_machine(DEFAULT_MACHINE)
THREE_D_SENTINEL = "__3d__"
MAX_SPEED_SAMPLES = 240
MAX_TORQUE_SAMPLES = 200
SURFACE_POINTS = 80
LOG = logging.getLogger(__name__)


def _parse_field_currents(raw: str, fallback: list[float]) -> list[float]:
    try:
        values = [float(val) for val in raw.replace(";", ",").split(",") if val.strip()]
        return values or fallback
    except Exception:
        return fallback


def _input_cell(label: str, component) -> html.Tr:
    return html.Tr([html.Th(label, style={"textAlign": "right", "paddingRight": "0.5rem"}), html.Td(component)])


def build_motor_tab() -> html.Div:
    machine_table = html.Table(
        [
            _input_cell("Vdc (V)", dcc.Input(id="motor-vdc", type="number", value=DEFAULT_MACHINE.Vdc, step=1)),
            _input_cell("Ls (H)", dcc.Input(id="motor-ls", type="number", value=DEFAULT_MACHINE.Ls, step=1e-4)),
            _input_cell("Lm (H)", dcc.Input(id="motor-lm", type="number", value=DEFAULT_MACHINE.Lm, step=1e-4)),
            _input_cell("Rs (Ω)", dcc.Input(id="motor-rs", type="number", value=DEFAULT_MACHINE.Rs, step=1e-4)),
            _input_cell("Rf (Ω)", dcc.Input(id="motor-rf", type="number", value=DEFAULT_MACHINE.Rf, step=0.1)),
            _input_cell("Poles", dcc.Input(id="motor-poles", type="number", value=DEFAULT_MACHINE.p, step=2, min=2)),
            _input_cell(
                "Electrical base ωₑ,b (rad/s)",
                dcc.Input(id="motor-wb", type="number", value=DEFAULT_MACHINE.wb, step=1),
            ),
            _input_cell(
                "Stator current limit |I₀|max (A pk)",
                dcc.Input(id="motor-i0max", type="number", value=DEFAULT_MACHINE.I0_max, step=1),
            ),
            _input_cell(
                "Field current limit |I_f|max (A dc)",
                dcc.Input(id="motor-ifmax", type="number", value=DEFAULT_MACHINE.If_max, step=0.1),
            ),
            _input_cell(
                "Field current cases (A, comma-separated)",
                dcc.Input(id="motor-field-currents", type="text", value="7.0,4.666,2.333"),
            ),
        ],
        style={"borderCollapse": "separate", "borderSpacing": "0.4rem 0.35rem", "width": "100%"},
    )

    plot_config = html.Table(
        [
            _input_cell(
                "Auto ranges",
                dcc.Checklist(
                    id="motor-auto-config",
                    options=[{"label": "Use defaults", "value": "auto"}],
                    value=["auto"],
                    style={"margin": "0"},
                ),
            ),
            _input_cell(
                "ωₘ max (rad/s)",
                dcc.Input(id="motor-speed-max", type="number", value=3 * DEFAULT_MACHINE.wmb, step=10),
            ),
            _input_cell(
                "Tₑ max (N·m)",
                dcc.Input(id="motor-torque-max", type="number", value=230, step=5),
            ),
            _input_cell(
                "Speed samples",
                dcc.Input(id="motor-speed-samples", type="number", min=40, max=240, step=10, value=140),
            ),
            _input_cell(
                "Torque samples",
                dcc.Input(id="motor-torque-samples", type="number", min=40, max=200, step=10, value=100),
            ),
        ],
        style={"borderCollapse": "separate", "borderSpacing": "0.4rem 0.35rem", "width": "100%"},
    )

    controls = html.Div(
        [
            html.Div([html.H4("Machine & inverter"), machine_table], className="motor-card"),
            html.Div([html.H4("Plot configuration"), plot_config], className="motor-card"),
        ],
        style={"display": "grid", "gridTemplateColumns": "repeat(auto-fit, minmax(360px, 1fr))", "gap": "0.75rem"},
    )

    action_row = html.Div(
        [
            html.Div(
                [
                    html.Button("Calculate", id="motor-calc-run", n_clicks=0, className="btn"),
                    html.Button("Download results", id="motor-download-btn", n_clicks=0, className="btn"),
                    dcc.Download(id="motor-download"),
                ],
                style={"display": "flex", "gap": "0.5rem", "flexWrap": "wrap"},
            ),
            html.Div(
                [
                    html.Label("Select I_f for contour plots"),
                    dcc.Dropdown(id="motor-selected-if", clearable=False),
                ],
                style={"minWidth": "240px"},
            ),
        ],
        style={"display": "flex", "justifyContent": "space-between", "alignItems": "center", "gap": "1rem"},
    )

    contour_graphs = html.Div(
        [
            dcc.Graph(id="motor-v0-graph"),
            dcc.Graph(id="motor-power-contour"),
            dcc.Graph(id="motor-modulation-contour"),
            dcc.Graph(id="motor-headroom"),
            dcc.Graph(id="motor-efficiency"),
            dcc.Graph(id="motor-loss-graph"),
            dcc.Graph(id="motor-loss-contour"),
            dcc.Markdown(
                id="motor-summary",
                mathjax=True,
                style={"backgroundColor": "#f8f9fa", "padding": "0.75rem", "borderRadius": "6px"},
            ),
        ],
        id="motor-contour-wrapper",
        style={"display": "grid", "gridTemplateColumns": "repeat(auto-fit, minmax(420px, 1fr))", "gap": "0.75rem"},
    )

    surface_graphs = html.Div(
        [
            dcc.Graph(id="motor-3d-power"),
            dcc.Graph(id="motor-3d-modulation"),
            dcc.Graph(id="motor-3d-headroom"),
            dcc.Graph(id="motor-3d-loss"),
            dcc.Graph(id="motor-3d-efficiency"),
        ],
        id="motor-surface-wrapper",
        style={"display": "none", "gridTemplateColumns": "repeat(auto-fit, minmax(420px, 1fr))", "gap": "0.75rem"},
    )

    return html.Div(
        [
            dcc.Store(id="motor-precomputed"),
            html.Div(id="motor-status-banner", style={"display": "flex", "alignItems": "center", "gap": "0.5rem"}),
            controls,
            action_row,
            dcc.Loading(contour_graphs, type="dot", id="motor-loading"),
            surface_graphs,
        ],
        className="motor-tab",
        style={"display": "flex", "flexDirection": "column", "gap": "0.75rem"},
    )


@callback(
    Output("motor-precomputed", "data"),
    Output("motor-selected-if", "options"),
    Output("motor-selected-if", "value"),
    Output("motor-status-banner", "children"),
    Output("app-header-status", "children"),
    Input("motor-calc-run", "n_clicks"),
    State("motor-vdc", "value"),
    State("motor-ls", "value"),
    State("motor-lm", "value"),
    State("motor-rs", "value"),
    State("motor-rf", "value"),
    State("motor-poles", "value"),
    State("motor-wb", "value"),
    State("motor-i0max", "value"),
    State("motor-ifmax", "value"),
    State("motor-field-currents", "value"),
    State("motor-auto-config", "value"),
    State("motor-speed-max", "value"),
    State("motor-torque-max", "value"),
    State("motor-speed-samples", "value"),
    State("motor-torque-samples", "value"),
    prevent_initial_call=True,
)
def compute_motor_maps(
    n_clicks,
    vdc,
    ls,
    lm,
    rs,
    rf,
    poles,
    wb,
    i0_max,
    if_max,
    field_current_text,
    auto_flags,
    speed_max,
    torque_max,
    speed_samples,
    torque_samples,
):
    del n_clicks  # unused for values
    start = time.perf_counter()

    params = MachineParams(
        Lm=lm,
        Ls=ls,
        Rs=rs,
        Rf=rf,
        p=int(poles),
        wb=wb,
        wmb=2.0 * wb / float(poles),
        Vdc=vdc,
        I0_max=i0_max,
        If_max=if_max,
    )
    bases = pu_bases_from_machine(params)
    speed_max = float(speed_max or 3 * params.wmb)
    torque_max = float(torque_max or bases.TB * 1.5)

    speed_samples = int(speed_samples or 0)
    torque_samples = int(torque_samples or 0)
    speed_samples = min(MAX_SPEED_SAMPLES, max(40, speed_samples))
    torque_samples = min(MAX_TORQUE_SAMPLES, max(40, torque_samples))

    auto_config = "auto" in (auto_flags or [])
    if auto_config:
        speed_max = 3 * params.wmb
        torque_max = max(params.I0_max * bases.TB / bases.IB, bases.TB * 1.5)

    omega_m = np.linspace(0.0, speed_max, speed_samples)
    torque_axis = np.linspace(0.0, torque_max, torque_samples)

    field_currents = _parse_field_currents(field_current_text or "", [if_max])
    options = ([{"label": "3D surfaces (all I_f)", "value": THREE_D_SENTINEL}] +
               [{"label": f"{val:.3f} A", "value": val} for val in field_currents])

    LOG.info(
        "Motor calc start: Vdc=%.1f, poles=%s, field_cases=%s, grid=%sx%s, auto=%s",
        vdc,
        poles,
        field_currents,
        speed_samples,
        torque_samples,
        auto_config,
    )

    voltage_sweep = inverter_voltage_vs_speed(params, bases, field_currents, omega_m)

    per_if_maps = {}
    for if_case in field_currents:
        power_map = mechanical_power_map(if_case, torque_axis, omega_m, params, bases)
        mod_map = modulation_index_map(if_case, torque_axis, omega_m, params, bases)
        loss_map = loss_surface(if_case, torque_axis, omega_m, params)
        eff_map = efficiency_map(if_case, torque_axis, omega_m, params, bases)

        per_if_maps[str(if_case)] = {
            "power": {k: v.tolist() for k, v in power_map.items()},
            "modulation": {k: v.tolist() for k, v in mod_map.items()},
            "loss": {k: v.tolist() for k, v in loss_map.items()},
            "efficiency": {k: v.tolist() for k, v in eff_map.items()},
        }

    payload = {
        "params": asdict(params),
        "bases": asdict(bases),
        "field_currents": field_currents,
        "omega_axis": omega_m.tolist(),
        "torque_axis": torque_axis.tolist(),
        "voltage_sweep": {
            "limit": voltage_sweep["limit"],
            "limit_pu": voltage_sweep["limit_pu"],
            "traces": [
                {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in trace.items()}
                for trace in voltage_sweep["traces"]
            ],
        },
        "per_if": per_if_maps,
    }

    payload_size = len(json.dumps(payload))
    duration = time.perf_counter() - start
    LOG.info(
        "Motor calc finished in %.2fs (%.1f kB, ω_max=%.1f rad/s, T_max=%.1f N·m)",
        duration,
        payload_size / 1024,
        speed_max,
        torque_max,
    )

    status = html.Span(
        [
            html.Strong("Motor calculator ready • "),
            f"ωₘ ≤ {speed_max:.1f} rad/s, Tₑ ≤ {torque_max:.1f} N·m across {len(field_currents)} field currents",
        ]
    )

    return payload, options, options[1]["value"], status, status


def _headroom_from_mod_map(mod_map: dict) -> np.ndarray:
    current_ratio = np.array(mod_map["current_ratio"], dtype=float)
    voltage_ratio = np.array(mod_map["voltage_ratio"], dtype=float)
    feasible_mask = np.array(mod_map["feasible_mask"], dtype=bool)

    headroom = np.minimum(
        np.divide(1.0, current_ratio, out=np.full_like(current_ratio, np.nan), where=current_ratio != 0),
        np.divide(1.0, voltage_ratio, out=np.full_like(voltage_ratio, np.nan), where=voltage_ratio != 0),
    )
    headroom = np.where(feasible_mask, headroom, np.nan)
    headroom = np.clip(headroom, 0, 5)
    return headroom


def _decimate_surface_arrays(x_grid: np.ndarray, y_grid: np.ndarray, z_grid: np.ndarray, metric: np.ndarray):
    row_step = max(1, math.ceil(metric.shape[0] / SURFACE_POINTS))
    col_step = max(1, math.ceil(metric.shape[1] / SURFACE_POINTS))
    return (
        x_grid[::row_step, ::col_step],
        y_grid[::row_step, ::col_step],
        z_grid[::row_step, ::col_step],
        metric[::row_step, ::col_step],
    )


@callback(
    Output("motor-v0-graph", "figure"),
    Output("motor-power-contour", "figure"),
    Output("motor-modulation-contour", "figure"),
    Output("motor-headroom", "figure"),
    Output("motor-efficiency", "figure"),
    Output("motor-loss-graph", "figure"),
    Output("motor-loss-contour", "figure"),
    Output("motor-summary", "children"),
    Output("motor-contour-wrapper", "style"),
    Output("motor-surface-wrapper", "style"),
    Output("motor-3d-power", "figure"),
    Output("motor-3d-modulation", "figure"),
    Output("motor-3d-headroom", "figure"),
    Output("motor-3d-loss", "figure"),
    Output("motor-3d-efficiency", "figure"),
    Input("motor-selected-if", "value"),
    Input("motor-precomputed", "data"),
)
def update_motor_plots(selected_if, data):
    if not data:
        placeholder = go.Figure()
        placeholder.update_layout(title="Run a calculation to populate plots")
        hidden = {"display": "none"}
        return (placeholder,) * 7 + ("", hidden, hidden, placeholder, placeholder, placeholder, placeholder, placeholder)

    contour_style = {"display": "grid", "gridTemplateColumns": "repeat(auto-fit, minmax(420px, 1fr))", "gap": "0.75rem"}
    surface_style = {"display": "none", "gridTemplateColumns": "repeat(auto-fit, minmax(420px, 1fr))", "gap": "0.75rem"}

    params = MachineParams(**data["params"])
    bases = PUBases(**data["bases"])
    omega_axis = np.array(data["omega_axis"], dtype=float)
    torque_axis = np.array(data["torque_axis"], dtype=float)
    field_currents = data["field_currents"]

    voltage_sweep = data["voltage_sweep"]
    v0_fig = go.Figure()
    for trace in voltage_sweep["traces"]:
        v0_fig.add_trace(
            go.Scatter(x=trace["omega_m"], y=trace["V0"], mode="lines", name=f"I_f = {trace['If']:.3f} A"),
        )
        if trace.get("cross_speed") is not None:
            v0_fig.add_trace(
                go.Scatter(
                    x=[trace["cross_speed"]], y=[voltage_sweep["limit"]], mode="markers", marker=dict(size=9), showlegend=False
                )
            )
    v0_fig.add_hline(
        y=voltage_sweep["limit"],
        line_dash="dash",
        line_color="#d62728",
        annotation_text=f"V₀ limit = Vdc/√3 = {voltage_sweep['limit']:.1f} V",
    )
    v0_fig.update_layout(title="Inverter V₀ vs Mechanical Speed", xaxis_title="ωₘ (rad/s)", yaxis_title="V₀ (Vₗₙ,peak)")

    # 3D mode branch
    if selected_if == THREE_D_SENTINEL:
        surface_style = contour_style
        contour_style = {"display": "none"}

        def surface_fig(metric_key: str, colorscale: str, color_title: str) -> go.Figure:
            fig = go.Figure()
            for if_case in field_currents:
                maps = data["per_if"][str(if_case)]
                omega_grid, torque_grid = np.meshgrid(omega_axis, torque_axis)
                y_grid = np.full_like(omega_grid, float(if_case))
                if metric_key == "power":
                    metric = np.array(maps["power"]["power_w"], dtype=float)
                elif metric_key == "modulation":
                    metric = np.array(maps["modulation"]["m0_masked"], dtype=float)
                elif metric_key == "headroom":
                    metric = _headroom_from_mod_map(maps["modulation"])
                elif metric_key == "loss":
                    metric = np.array(maps["loss"]["loss_kw"], dtype=float)
                else:
                    metric = np.array(maps["efficiency"]["efficiency"], dtype=float)

                finite_metric = metric[np.isfinite(metric)]
                cmin = float(finite_metric.min()) if finite_metric.size else 0.0
                cmax = float(finite_metric.max()) if finite_metric.size else 1.0

                omega_ds, y_ds, torque_ds, metric_ds = _decimate_surface_arrays(
                    omega_grid, y_grid, torque_grid, metric
                )

                fig.add_trace(
                    go.Surface(
                        x=omega_ds,
                        y=y_ds,
                        z=torque_ds,
                        surfacecolor=metric_ds,
                        colorscale=colorscale,
                        cmin=cmin,
                        cmax=cmax,
                        showscale=if_case == field_currents[-1],
                        colorbar=dict(title=color_title) if if_case == field_currents[-1] else None,
                        opacity=0.9,
                        name=f"I_f={if_case:.3f} A",
                    )
                )
            fig.update_layout(
                scene=dict(
                    xaxis_title="ωₘ (rad/s)",
                    yaxis_title="I_f (A)",
                    zaxis_title="Tₑ (N·m)",
                    bgcolor="white",
                ),
                coloraxis_colorbar=dict(title=color_title),
                margin=dict(l=0, r=0, b=0, t=30),
            )
            return fig

        return (
            v0_fig,
            go.Figure(),
            go.Figure(),
            go.Figure(),
            go.Figure(),
            go.Figure(),
            go.Figure(),
            "",
            contour_style,
            surface_style,
            surface_fig("power", "Viridis", "Pₘ"),
            surface_fig("modulation", "Plasma", "M₀"),
            surface_fig("headroom", "Cividis", "Headroom"),
            surface_fig("loss", "Magma", "P_loss"),
            surface_fig("efficiency", "Blues", "η"),
        )

    maps = data["per_if"].get(str(selected_if))
    if maps is None:
        return no_update

    power_map = {k: np.array(v, dtype=float) for k, v in maps["power"].items()}
    mod_map = {k: np.array(v, dtype=float) for k, v in maps["modulation"].items()}
    loss_map = {k: np.array(v, dtype=float) for k, v in maps["loss"].items()}
    eff_map = {k: np.array(v, dtype=float) for k, v in maps["efficiency"].items()}

    power_fig = go.Figure(
        go.Contour(
            x=power_map["omega_m"][0, :],
            y=power_map["torque"][:, 0],
            z=power_map["power_w"] / 1e3,
            colorscale="Viridis",
            colorbar=dict(title="Pₘ (kW)"),
            contours=dict(showlabels=True),
        )
    )
    power_fig.add_trace(
        go.Contour(
            x=power_map["omega_m"][0, :],
            y=power_map["torque"][:, 0],
            z=power_map["voltage_mask"].astype(float),
            contours=dict(coloring="lines", showlines=True, start=1, end=1, size=1),
            showscale=False,
            line=dict(color="#d62728"),
            name="V₀ = Vdc/√3",
        )
    )
    power_fig.add_trace(
        go.Contour(
            x=power_map["omega_m"][0, :],
            y=power_map["torque"][:, 0],
            z=power_map["current_mask"].astype(float),
            contours=dict(coloring="lines", showlines=True, start=1, end=1, size=1),
            showscale=False,
            line=dict(color="#7f7f7f", dash="dash"),
            name="|I| = I₀,max",
        )
    )
    power_fig.update_layout(
        title=f"Mechanical Power Contours | I_f = {float(selected_if):.3f} A",
        xaxis_title="ωₘ (rad/s)",
        yaxis_title="Tₑ (N·m)",
    )

    mod_fig = go.Figure(
        go.Contour(
            x=mod_map["omega_m"][0, :],
            y=mod_map["torque"][:, 0],
            z=mod_map["m0_masked"],
            colorscale="Plasma",
            colorbar=dict(title="M₀"),
            contours=dict(showlabels=True),
        )
    )
    mod_fig.add_trace(
        go.Contour(
            x=mod_map["omega_m"][0, :],
            y=mod_map["torque"][:, 0],
            z=mod_map["m0"],
            contours=dict(coloring="lines", showlines=True, start=2 / math.sqrt(3.0), end=2 / math.sqrt(3.0), size=1),
            showscale=False,
            line=dict(color="#d62728"),
            name="Linear modulation limit",
        )
    )
    mod_fig.update_layout(
        title=f"Modulation Index Contours | I_f = {float(selected_if):.3f} A",
        xaxis_title="ωₘ (rad/s)",
        yaxis_title="Tₑ (N·m)",
    )

    headroom = _headroom_from_mod_map(maps["modulation"])
    headroom_max = float(np.nanmax(headroom)) if np.isfinite(headroom).any() else 1.0
    headroom_fig = go.Figure(
        go.Contour(
            x=mod_map["omega_m"][0, :],
            y=mod_map["torque"][:, 0],
            z=headroom,
            colorscale="Cividis",
            colorbar=dict(title="Headroom (× limit)"),
            contours=dict(showlabels=True),
            zmin=0,
            zmax=headroom_max,
        )
    )
    headroom_fig.add_trace(
        go.Contour(
            x=mod_map["omega_m"][0, :],
            y=mod_map["torque"][:, 0],
            z=mod_map["voltage_mask"].astype(float),
            contours=dict(coloring="lines", showlines=True, start=1, end=1, size=1),
            showscale=False,
            line=dict(color="#d62728"),
            name="Voltage limit",
        )
    )
    headroom_fig.add_trace(
        go.Contour(
            x=mod_map["omega_m"][0, :],
            y=mod_map["torque"][:, 0],
            z=mod_map["current_mask"].astype(float),
            contours=dict(coloring="lines", showlines=True, start=1, end=1, size=1),
            showscale=False,
            line=dict(color="#7f7f7f", dash="dash"),
            name="Current limit",
        )
    )
    headroom_fig.update_layout(
        title=f"Voltage/Current Headroom | I_f = {float(selected_if):.3f} A",
        xaxis_title="ωₘ (rad/s)",
        yaxis_title="Tₑ (N·m)",
    )

    eff_fig = go.Figure(
        go.Contour(
            x=eff_map["omega_m"][0, :],
            y=eff_map["torque"][:, 0],
            z=eff_map["efficiency"],
            colorscale="Blues",
            colorbar=dict(title="η"),
            contours=dict(showlabels=True),
        )
    )
    eff_fig.update_layout(
        title=f"Efficiency Contours | I_f = {float(selected_if):.3f} A",
        xaxis_title="ωₘ (rad/s)",
        yaxis_title="Tₑ (N·m)",
    )

    losses = loss_curves(params)
    loss_fig = go.Figure()
    loss_fig.add_trace(go.Scatter(x=losses["If"], y=losses["Pf"], name="P_f"))
    loss_fig.add_trace(go.Scatter(x=losses["E0"], y=losses["Pcore"], name="P_core"))
    loss_fig.add_trace(go.Scatter(x=losses["I0"], y=losses["Ps_total"], name="P_s + P_inv"))
    loss_fig.update_layout(title="Loss building blocks", xaxis_title="Input sweep", yaxis_title="Loss (W)")

    loss_contour = go.Figure(
        go.Contour(
            x=loss_map["omega_m"][0, :],
            y=loss_map["torque"][:, 0],
            z=loss_map["loss_kw"],
            colorscale="Magma",
            colorbar=dict(title="P_loss (kW)"),
            contours=dict(showlabels=True),
        )
    )
    loss_contour.update_layout(
        title=f"Total Loss Contours | I_f = {float(selected_if):.3f} A",
        xaxis_title="ωₘ (rad/s)",
        yaxis_title="Tₑ (N·m)",
    )

    base_summary = (
        f"**Per-unit bases**\n"
        f"* $V_B = V_{{dc}}/\\sqrt{{3}} = {bases.VB:.1f}$ V (LN peak)\n"
        f"* $I_B = {bases.IB:.1f}$ A pk, $Z_B = {bases.ZB:.3f}$ Ω, $L_B = {bases.LB:.6f}$ H\n"
        f"* $\\omega_{{m,B}} = {bases.wmB:.1f}$ rad/s, $T_B = {bases.TB:.1f}$ N·m\n\n"
        f"**Selected case**: $I_f = {float(selected_if):.3f}$ A, range $\\omega_m \\le {omega_axis.max():.1f}$ rad/s, "
        f"$T_e \\le {torque_axis.max():.1f}$ N·m."
    )

    return (
        v0_fig,
        power_fig,
        mod_fig,
        headroom_fig,
        eff_fig,
        loss_fig,
        loss_contour,
        base_summary,
        contour_style,
        surface_style,
        go.Figure(),
        go.Figure(),
        go.Figure(),
        go.Figure(),
        go.Figure(),
    )


@callback(
    Output("motor-download", "data"),
    Input("motor-download-btn", "n_clicks"),
    State("motor-precomputed", "data"),
    prevent_initial_call=True,
)
def download_results(n_clicks, data):
    del n_clicks
    if not data:
        return no_update

    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("motor_results.json", json.dumps(data))
    buffer.seek(0)
    return dcc.send_bytes(buffer.read, filename="motor_results.zip")
