"""Tests for the access-token gate on the Dash server."""

import importlib


def test_static_assets_bypass_token_guard(monkeypatch):
    import app as app_module

    monkeypatch.setenv("ACCESS_TOKEN", "secret-token")
    reloaded = importlib.reload(app_module)

    client = reloaded.app.server.test_client()

    # Static assets should not be blocked by the access token guard.
    static_response = client.get("/_dash-component-suites/fake.js")
    assert static_response.status_code != 403

    # The main page still requires the token in the query string.
    assert client.get("/").status_code == 403
    assert client.get("/?token=secret-token").status_code == 200

    # Reload the module with the token removed to avoid side-effects on other tests.
    monkeypatch.delenv("ACCESS_TOKEN", raising=False)
    importlib.reload(app_module)
