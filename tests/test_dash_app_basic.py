import shutil

import pytest
from selenium.webdriver import ActionChains

from app import app, APP_TITLE


pytestmark = pytest.mark.skipif(
    shutil.which("chromedriver") is None,
    reason="chromedriver is required for dash[testing] browser fixtures",
)


def test_app_layout_has_tabs(dash_duo):
    dash_duo.start_server(app)
    dash_duo.wait_for_text_to_equal("h1", APP_TITLE, timeout=10)
    tabs = dash_duo.find_elements("div.tab")
    assert len(tabs) >= 4


@pytest.mark.parametrize("offset", [60])
def test_theta_slider_updates_summary(dash_duo, offset):
    dash_duo.start_server(app)
    dash_duo.wait_for_text_to_equal("h1", APP_TITLE, timeout=10)
    summary = dash_duo.find_element("#frames-summary").text
    handle = dash_duo.wait_for_element("#frames-theta .rc-slider-handle")
    ActionChains(dash_duo.driver).drag_and_drop_by_offset(handle, offset, 0).perform()
    dash_duo.wait_until(lambda: dash_duo.find_element("#frames-summary").text != summary, timeout=4)
