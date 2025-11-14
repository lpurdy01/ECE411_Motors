# Repository Guidance for AI Contributors

Welcome! This project powers an interactive Dash application used to build
intuition for ECE 411 synchronous machine concepts. When making changes inside
this repository, please follow the conventions below so that future
contributors (human or AI) can work together smoothly.

## General Engineering Practices

1. **Preserve interactivity.** The Dash UI is the primary teaching tool. When
   editing UI code, prefer descriptive controls, responsive plots, and rich
   textual explanations. Favor `dcc.Markdown` for math-heavy narratives so we
   can embed LaTeX snippets.
2. **Keep math reusable.** All machine calculations belong in
   `motor_models.py`. This module must remain free of Dash imports and rely on
   NumPy for vectorized math so tests stay fast.
3. **Document formulas.** Whenever you implement a new equation, add an inline
   comment or docstring that shows both the symbolic expression (LaTeX is ok)
   and an example of how variables plug into it.
4. **Tests first (or at least together).** New math helpers or callbacks should
   ship with pytest coverage in `tests/`. If a change alters existing behavior,
   update the relevant test expectations.
5. **Lint-by-principle.** We do not run a dedicated linter yet, but please keep
   imports sorted, favor explicit names, and avoid wildcard imports.

## Dash / UI Specifics

* Prefer layout components from `dash.dcc` and `dash.html`; avoid raw HTML
  strings.
* For textual summaries, use `dcc.Markdown` with `mathjax=True` so equations
  render properly.
* Keep callbacks pure: no hidden global state beyond module-level constants
  such as sample grids.

## Testing & Tooling

* Run `pytest` locally before committing. Dash smoke tests skip automatically
  when `chromedriver` is unavailable.
* GitHub Actions in `.github/workflows/ci.yml` must stay green—update it if new
  dependencies are introduced.

Thanks for helping maintain a clear, educational experience!
