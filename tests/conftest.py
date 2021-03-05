import os

import pytest
# import chromedriver_binary  # Overwrite "chromedriver" installed by the package with one of a version matching Chrome installed in the machine
from selenium import webdriver
from selenium.common.exceptions import StaleElementReferenceException
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.desired_capabilities import DesiredCapabilities
import sys
import datetime
import time
import ctypes.util
from typing import Union
# from pyvirtualdisplay import Display

# scopes (broader to narrower): session, module, class, function (default)
from biobarcoding.rest.main import create_app
from pathlib import Path


@pytest.fixture(scope="module")
def testful():
    home = str(Path.home())
    cfg = dict(DB_CONNECTION_STRING="postgresql://postgres:postgres@localhost:5432/",
               CACHE_FILE_LOCATION=f"{home}/.cache/bcs-backend/cache",
               REDIS_HOST_FILESYSTEM_DIR=f"{home}/.cache/bcs-backend/sessions",
               REDIS_HOST="filesystem:local_session",
               TESTING="True",
               SELF_SCHEMA="")
    app = create_app(True, cfg)
    testing_client = app.test_client()
    ctx = app.app_context()
    ctx.push()

    yield testing_client

    ctx.pop()
    # TODO - ISSUE: "after_a_request" function seems to be not called. WHY? SOLVE


@pytest.fixture(scope="session")
def chrome_driver():
    """
    Prepare Chrome driver, visible or headless
    :return:
    """
    options = Options()
    options.add_argument("--disable-extensions")
    options.add_experimental_option("useAutomationExtension", False)
    options.binary_location = "/opt/google/chrome/chrome"
    # options.headless = True
    driver = webdriver.Chrome(options=options)

    yield driver # Execute tests

    driver.quit()


@pytest.fixture(scope="session")
def firefox_driver():
    """
    Prepare Firefox driver, visible or headless
    :return:
    """
    caps = DesiredCapabilities.FIREFOX
    caps["wires"] = True
    options = webdriver.FirefoxOptions()
    # options.headless = True
    # display = Display(visible=0, size=(1024, 769))
    # display.start()

    driver = webdriver.Firefox(options=options, capabilities=caps)

    yield driver # Execute tests

    driver.quit()


@pytest.fixture(scope="session")
def base_data():
    # TODO Load base data (using dbUnit or other)
    return


# ----------------------------------------------------------------------------------------------------------------------
# UTILS
# ----------------------------------------------------------------------------------------------------------------------

def find_elements_through_shadow_roots(dr, path: str):
    """
    Obtain all the elements matching the "enhanced" XPath passed in "path"

    "path" is an XPath elaborated following these steps:

    * In Chrome Developer Tools, take the XPath of the desired element, right click menu over the desired element, then "Copy -> Copy full XPath"
    * If there are "shadow DOM" elements in the middle, "//" in the Xpath, substitute with "/#shadow-root/"
    * The part after a shadow DOM element can use a CSS selector as obtained by "querySelectorAll", and then an XPath predicate between square brackets "[...]"
    * ADVISE: in the original "Copy full XPath", avoid indexed elements, try to find a property not depending on the exact order

    :param dr: WebDriver instance
    :param path: A string following the previous description
    :return:
    """
    def find_elements_in_shadowroot(dr, element, child_name):
        s = """
    // document.evaluate('app-header-layout/app-header/app-toolbar/header-menu', element, null, XPathResult.ANY_TYPE, null ).iterateNext()    
    // document.evaluate('/html/body/vaadin-context-menu-overlay', document, null, XPathResult.ANY_TYPE, null ).iterateNext()
    return arguments[0].shadowRoot.querySelectorAll(arguments[1]); 
        """
        return dr.execute_script(s, element, child_name)

    def prepare_parts():
        import re
        parts = []
        for i, xpath in enumerate(re.split("#shadow-root", path)):
            if i > 0:
                # Split in two, and update "xpath" with "the rest"
                lst = xpath[1:].split("/", 1)
                if len(lst) > 1:
                    xpath = lst[1]
                else:
                    xpath = ""
                # Split if "[" present (needs XPath)
                if "[" in lst[0]:
                    pos = lst[0].find("[")
                    part = lst[0][0:pos]
                    rest = lst[0][pos:]
                else:
                    part = lst[0]
                    rest = None
                parts.append(("shadow-root", part))
                # Partial XPath
                if rest:
                    parts.append(("xpath", f"*{rest}"))

            if xpath != "" and xpath != "/":
                if xpath.endswith("/"):
                    xpath = xpath[:-1]

                parts.append(("xpath", xpath))
        return parts

    def evaluate_part(elem, part):
        if part[0] == "xpath":
            if elem:
                return elem.find_elements_by_xpath(part[1])
            else:
                return dr.find_elements_by_xpath(part[1])
        else:
            return find_elements_in_shadowroot(dr, elem, part)

    def process_level(starting_element, parts, level):
        tmp = evaluate_part(starting_element, parts[level])
        print(f"Level: {level}, Part: {parts[level]}, # candidate elements: {len(tmp)}")
        if level == len(parts) - 1:
            print(f"Last level, length {len(tmp)}")
            return tmp
        else:
            ret = []
            for el in tmp:
                ret.extend(process_level(el, parts, level + 1))
            print(f"Finished recursive level {level}, length {len(ret)}")
            return ret

    return process_level(None, prepare_parts(), 0)


def find_element_through_shadow_roots(dr, path: str):
    r = find_elements_through_shadow_roots(dr, path)
    if len(r) != 1:
        raise Exception(f"Se encontraron ({len(r)} que respondieron a la ruta '{path}'. Debe ser un único elemento.")
    return r[0]


class CheckPresenceOfElementAtPath(object):
    """ An expectation for checking that an element is present on the DOM
    of a page. This does not necessarily mean that the element is visible.
    locator - used to find the element
    returns the WebElement once it is located
    """
    def __init__(self, path):
        self._path = path

    def __call__(self, driver):
        return find_element_through_shadow_roots(driver, self._path)
