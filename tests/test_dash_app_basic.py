from dash import dcc

from app import app


def test_app_layout_has_expected_tabs():
    layout = app.layout
    tabs_component = next(child for child in layout.children if isinstance(child, dcc.Tabs))
    assert len(tabs_component.children) == 6
    labels = [tab.label for tab in tabs_component.children]
    assert labels == [
        "Frame Playground",
        "Round-Rotor",
        "Field Weakening",
        "Salient Torque",
        "Motor Calculator",
        "Cheat Sheet",
    ]


def test_index_page_served_successfully():
    client = app.server.test_client()
    response = client.get("/")
    assert response.status_code == 200
    assert "DashRenderer" in response.get_data(as_text=True)
