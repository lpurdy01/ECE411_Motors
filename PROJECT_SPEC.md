# Project Summary – ECE 411 Motor Visualization Dash App

## 1. Purpose

Build a **Python Dash web app** that helps visualize and understand:

* Clarke & Park transforms (abc → αβ → qd0) under **Lipo–Novotny** conventions.
* Synchronous machine phasors and torque (round-rotor and salient-pole).
* Field-weakening, voltage/current limits, and torque–speed behavior.

The app is an **interactive learning tool** aimed at ECE 411-style homework and exam concepts (HW3 & HW4).

Target runtime: **Python 3.10+** in a GitHub Codespace or local venv.
Primary stack: **Dash + Plotly + NumPy + pytest + GitHub Actions**.

---

## 2. High-Level Features (Definition of “Complete Implementation”)

A *complete implementation* should include:

1. **Modular core math library** (`motor_models.py`)

   * Implements all transforms and machine equations in pure Python/NumPy.
   * No Dash/UI dependencies in this layer.
2. **Dash app with multiple tabs** (`app.py`)

   * Tab 1: abc ↔ αβ ↔ qd0 “Frame Playground”
   * Tab 2: Round-rotor phasor & torque
   * Tab 3: Field-weakening & torque–speed
   * Tab 4: Salient-pole torque vs δ (field + reluctance)
3. **Unit tests** (`tests/`) for:

   * Transforms (Clarke, Park).
   * Torque and power equations.
   * Field-weakening formulas.
4. **Integration tests** for:

   * Dash app startup and layout (tabs present, key components exist).
   * Basic callback behavior (changing inputs affects outputs).
5. **GitHub Actions CI** (`.github/workflows/ci.yml`)

   * Install dependencies.
   * Run `pytest`.
   * Optionally run basic Dash integration tests in headless mode.

If those are working and tests are green in CI, call the initial version “done”.

---

## 3. Directory Structure

Proposed repo layout:

```text
motor-viz/
  app.py                     # Dash app entry point
  motor_models.py            # Core math & model functions (pure Python)
  callbacks.py               # Optional: split Dash callbacks out of app.py
  requirements.txt           # Runtime dependencies
  dev-requirements.txt       # Dev/test deps (pytest, dash[testing], etc.)
  README.md                  # High-level overview
  PROJECT_SPEC.md            # This spec (optional)

  tests/
    test_transforms.py       # Unit tests for abc→αβ→qd0 and back
    test_round_rotor.py      # Unit tests for round-rotor equations
    test_field_weakening.py  # Unit tests for FW formulas
    test_salient_pole.py     # Unit tests for field/reluctance torque
    test_dash_app_basic.py   # Simple integration tests for Dash app

  .github/
    workflows/
      ci.yml                 # GitHub Actions workflow

  .devcontainer/             # (Optional) Codespaces config
    devcontainer.json
    Dockerfile
```

---

## 4. Core Math API (`motor_models.py`)

Implement a clean, testable API.

### 4.1 Clarke & Park transforms

```python
import numpy as np

def abc_to_alphabeta(va, vb, vc):
    """Return (alpha, beta, v0) using Lipo–Novotny Clarke transform."""
    alpha = (2/3) * (va - 0.5*vb - 0.5*vc)
    beta  = (2/3) * (np.sqrt(3)/2) * (vc - vb)  # note (Vc - Vb)
    v0    = (va + vb + vc) / 3
    return alpha, beta, v0

def alphabeta_to_qd(alpha, beta, theta_e):
    """Return (vq, vd) using Lipo–Novotny Park transform."""
    c = np.cos(theta_e)
    s = np.sin(theta_e)
    vq =  alpha * c - beta * s
    vd =  alpha * s + beta * c
    return vq, vd

def abc_to_qd0(va, vb, vc, theta_e):
    alpha, beta, v0 = abc_to_alphabeta(va, vb, vc)
    vq, vd = alphabeta_to_qd(alpha, beta, theta_e)
    return vq, vd, v0
```

### 4.2 Round-rotor synchronous machine

Assume per-unit:

```python
def round_rotor_voltage_current(v0, rs, ls, w_e, e0, delta):
    """
    Compute Vq,Vd, Iq,Id, I0, P, Te for a round-rotor SM.
    Lipo–Novotny dq, per-unit, steady-state.
    Inputs:
        v0   : terminal magnitude (pu)
        rs   : stator resistance (pu)
        ls   : synchronous inductance (pu)
        w_e  : electrical speed (pu)
        e0   : internal emf magnitude (pu)
        delta: load angle between E and V (rad)
    """
    # V phasor (qd)
    vq = v0 * np.cos(delta)
    vd = v0 * np.sin(delta)

    eq = e0  # E_d = 0 for round-rotor
    ed = 0.0

    # Solve linear system:
    # vq = rs*iq + w_e*ls*id + eq
    # vd = rs*id - w_e*ls*iq + ed
    # ...
    # Return Iq,Id, magnitude, P, Te, etc.
    ...
```

### 4.3 Salient-pole torque vs δ

```python
def salient_pole_torque_curves(ld, lq, e0, v0, w_e, delta_array):
    """
    Compute field torque, reluctance torque, and total torque vs delta.
    Uses:
      T_field = -(E0*V0)/(w_e*Ld) * sin(delta)
      T_rel   = -(V0**2*(Ld-Lq))/(2*w_e*Ld*Lq) * sin(2*delta)
    Returns arrays: (Te_field, Te_rel, Te_total)
    """
    t_field = -(e0 * v0)/(w_e * ld) * np.sin(delta_array)
    t_rel   = -(v0**2 * (ld - lq))/(2 * w_e * ld * lq) * np.sin(2 * delta_array)
    t_total = t_field + t_rel
    return t_field, t_rel, t_total
```

### 4.4 Field-weakening helpers

```python
def field_weakening_flux(v0, i0, ls, w_e):
    """
    R_s = 0 case:
      V0^2 = (w_e*Lambda0)^2 - (w_e*Ls*I0)^2
      Lambda0 = sqrt((V0/w_e)**2 + (Ls*I0)**2)
    """
    return np.sqrt((v0 / w_e)**2 + (ls * i0)**2)

def constant_power_torque(v0, i0, w_e):
    """Te_max and P_max in the constant-power region."""
    p_max = v0 * i0
    te_max = p_max / w_e
    return te_max, p_max
```

---

## 5. Dash App Requirements (`app.py`)

Use **Dash with tabs**. Structure:

* `dcc.Tabs` with 4 `dcc.Tab` children:

  1. “Frames”
  2. “Round-Rotor”
  3. “Field Weakening”
  4. “Salient Torque”

Each tab:

* A **control panel** (`html.Div` with `dcc.Slider`, `dcc.Input`, `dcc.Dropdown`).
* One or more `dcc.Graph` components with Plotly visualizations.
* A small text area (`html.Div`) for explaining what is currently shown.

Minimal requirement per tab:

### Tab 1 – Frames

* Sliders for:

  * `Vpk`, `f`, θ₀, harmonic order h, harmonic ratio.
  * θₑ (manual) or toggle for θₑ = ωₑt.
* Graphs:

  * Time-domain Va, Vb, Vc.
  * αβ vs time.
  * Phasor-style 2D plot of vαβ.
  * q/d vs time, given θₑ.

### Tab 2 – Round-Rotor

* Inputs: `Rs`, `Ls`, `V0`, `E0`, `w_e`, `delta`.
* Graph:

  * qd-plane phasor diagram (V, E, I, jωLsI, RsI).
* Text: P, Te, I0, power factor.

### Tab 3 – Field Weakening

* Inputs: `V0_max`, `I0_max`, `Ls`, `w_e` slider.
* Graphs:

  * Te vs ω_e.
  * P vs ω_e.
  * Λ0 vs ω_e (computed by helper).
* Show constant-torque and constant-power regions clearly.

### Tab 4 – Salient Torque

* Inputs: `Ld`, `Lq`, `E0`, `V0`, `w_e`.
* Graph:

  * T_field(δ), T_rel(δ), T_total(δ) vs δ for δ in [−π, π].
* Text: `T_field_max`, `T_rel_max`, saliency ratio.

---

## 6. Testing Strategy

Use **pytest** for all tests. `dev-requirements.txt` should include:

```text
pytest
dash[testing]
```

### 6.1 Unit Tests

#### `test_transforms.py`

* Test that a balanced abc set with known phase shift gives expected α, β, q, d for simple cases.

Examples:

* Case: (Va = cos(0), Vb = cos(-2π/3), Vc = cos(-4π/3)) at θₑ = 0:

  * Expect α ≈ 1, β ≈ 0.
  * Expect q ≈ 1, d ≈ 0.

* Case with θₑ = π/2:

  * Expect q ≈ 0, d ≈ 1 (with sign per convention).

#### `test_round_rotor.py`

* Give simple parameters (e.g., Rs=0, Ls=1, V0=1, E0 values).
* For δ = 0, check that I_d ≈ 0 and torque matches analytical formula.

#### `test_field_weakening.py`

* Check `field_weakening_flux` against hand calculation.
* Check that Te_max * ω_e = V0 * I0 (within tolerance).

#### `test_salient_pole.py`

* For chosen Ld, Lq, E0, V0, w_e:

  * Check that `T_field_max` ~ (E0*V0)/(w_e*Ld).
  * Check that `T_rel_max` ~ V0^2*(Ld-Lq)/(2*w_e*Ld*Lq).

### 6.2 Dash Integration Tests (`test_dash_app_basic.py`)

Use `dash[testing]` and its `dash_duo` fixture.

Basic smoke tests:

```python
import pytest
from app import app

def test_app_loads(dash_duo):
    dash_duo.start_server(app)
    dash_duo.wait_for_text_to_equal("h1", "ECE 411 Motor Visualization", timeout=10)
    tabs = dash_duo.find_elements("div .tab")
    assert len(tabs) >= 4


def test_frames_tab_updates(dash_duo):
    dash_duo.start_server(app)
    # Example: find the slider for theta_e and move it
    # Then assert that one of the graphs has updated (by checking plotly JSON length)
    ...
```

Minimum expectations:

* App starts without exceptions.
* Tabs exist and are clickable.
* Changing a control (e.g., slider) triggers at least one graph update.

---

## 7. GitHub Actions CI (`.github/workflows/ci.yml`)

A simple CI workflow:

```yaml
name: CI

on:
  push:
  pull_request:

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - name: Set up Python
        uses: actions/setup-python@v5
        with:
          python-version: '3.11'

      - name: Install dependencies
        run: |
          python -m pip install --upgrade pip
          pip install -r requirements.txt
          pip install -r dev-requirements.txt

      - name: Run tests
        run: pytest
```

If `dash[testing]` proves flaky in GitHub Actions, you can:

* Keep all **math unit tests** in CI.
* Mark Dash tests as `@pytest.mark.slow` and optionally skip in CI.

---

## 8. Codespaces / Devcontainer (Optional)

Provide a basic `.devcontainer/devcontainer.json`:

```json
{
  "name": "motor-viz",
  "image": "mcr.microsoft.com/devcontainers/python:3.11",
  "postCreateCommand": "pip install -r requirements.txt && pip install -r dev-requirements.txt",
  "forwardPorts": [8050]
}
```

This lets you hit `http://localhost:8050` in the Codespaces forwarded port to use the Dash app.

---

If you want, I can next generate a starter `app.py` and `motor_models.py` skeleton aligned with this spec so Copilot has a strong pattern to follow.
