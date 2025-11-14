# ECE411 Motors Visualization Toolkit

This repository now hosts an interactive **Dash web application** that turns key
ECE 411 synchronous machine equations into visual, hands-on demonstrations.
Use it to explore Clarke/Park transforms, round-rotor phasors, field-weakening,
and salient-pole torque contributions.

## Getting started

1. Create and activate a Python 3.10+ virtual environment.
2. Install dependencies:

   ```bash
   pip install -r requirements.txt
   # For development and testing
   pip install -r dev-requirements.txt
   ```

3. Launch the Dash app:

   ```bash
   python app.py
   ```

   The server defaults to `http://127.0.0.1:8050/`.

## Running tests

All unit and integration tests use `pytest` and `dash[testing]`:

```bash
pytest
```

The GitHub Actions workflow in `.github/workflows/ci.yml` mirrors this command.

## Project specification

The full project statement, including feature requirements and testing plan, is
captured in [`PROJECT_SPEC.md`](PROJECT_SPEC.md). It reflects the learning goals
for ECE 411 homework sets covering synchronous machine analysis.

## Repository layout

```
app.py                # Dash entry point with all layouts and callbacks
motor_models.py       # Pure-Python math helpers (transforms, torque, limits)
requirements.txt      # Runtime dependencies
dev-requirements.txt  # Test/development dependencies
PROJECT_SPEC.md       # High-level design and requirements
.github/workflows/    # Continuous integration configuration
tests/                # Unit + integration tests
```

Historical MATLAB scripts and reference material from the original coursework
are still included for reference under the existing directories.
