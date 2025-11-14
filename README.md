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

## Usage

Visit the locally running app at the address above, or the deployed service at
`https://<service>.onrender.com` (append `?token=YOURTOKEN` if you've set
`ACCESS_TOKEN`).

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

## Deployment via Render

This app is deployed automatically on commit to `main` branch.

1. On Render.com, create a new **Web Service**:
   - Link your GitHub repo.
   - Select branch `main`.
   - Build Command: `pip install -r requirements.txt`
   - Start Command: `gunicorn app:server`

2. In Render service settings → Environment Variables:
   - Add `ACCESS_TOKEN` (optional) with your secret value if you want token access.

3. Ensure auto-deploy is set to **On Commit** (default) so that every merge to `main` triggers a new deploy.  
   See Render docs on auto-deploy. :contentReference[oaicite:1]{index=1}

4. You will get a service URL like `https://<name>.onrender.com`.  
   Share this URL (plus `?token=…` if using token) with class members.
