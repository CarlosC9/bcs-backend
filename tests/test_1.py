"""
With "pytest", tests are functions starting with "test_". Their parameters are injected

<FUNCTIONS TAKEN FROM "SGI TESTING" PROJECT AS EXAMPLES>

"""
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions
from selenium.webdriver.support.wait import WebDriverWait

from .conftest import firefox_driver, find_element_through_shadow_roots, CheckPresenceOfElementAtPath

delay = 8

base_url = "http://rnebot-pc:8090"


def test_open_login_page(firefox_driver, time_set):
    """ LOAD INITIAL BCS PAGE AND CHECK SOMETHING. IF IT WORKS, THE SYSTEM IS RUNNING """
    driver = firefox_driver
    # set_system_date("2019-10-11T15:30:45")
    driver.get(base_url)
    WebDriverWait(driver, delay).until(expected_conditions.presence_of_element_located((By.XPATH, "/html/body/login-view/div/div/h1")))
    elem = driver.find_element_by_xpath("/html/body/login-view/div/div/h1")
    assert elem.text == ""


def load_and_login(d):
    d.get(base_url)
    usuario_label_xpath = '/html/body/login-view/div/iron-form/form[1]/div/vaadin-text-field/div/label'
    WebDriverWait(d, delay).until(expected_conditions.presence_of_element_located((By.XPATH, usuario_label_xpath)))
    elem = d.find_element_by_xpath(usuario_label_xpath)
    assert elem.text == "Usuario"
    # Login
    user = d.find_element_by_xpath('/html/body/login-view/div/iron-form/form[1]/div/vaadin-text-field/div/div[1]/input')
    user.send_keys("admin@example.com")
    passw = d.find_element_by_xpath('/html/body/login-view/div/iron-form/form[1]/div/vaadin-password-field/div/div[1]/input')
    passw.send_keys("admin")
    button = d.find_element_by_xpath('/html/body/login-view/div//form/div/vaadin-button/button')
    button.click()
    WebDriverWait(d, delay).until(expected_conditions.presence_of_element_located((By.XPATH, "/html/body/main-view")))


def logout(d):
    button = find_element_through_shadow_roots(d, "/html/body/main-view/#shadow-root/app-drawer-layout/app-header-layout/app-header/app-toolbar/header-menu/#shadow-root/vaadin-context-menu/vaadin-button/#shadow-root/button")
    button.click()
    p_link = "/html/body/vaadin-context-menu-overlay/#shadow-root/div#overlay/div/#shadow-root/div/vaadin-list-box/vaadin-item[contains(@class,'footer')]/a"
    WebDriverWait(d, delay).until(CheckPresenceOfElementAtPath(p_link))
    link = find_element_through_shadow_roots(d, p_link)
    link.click()
    usuario_label_xpath = '/html/body/login-view/div/iron-form/form[1]/div/vaadin-text-field/div/label'
    WebDriverWait(d, delay).until(expected_conditions.presence_of_element_located((By.XPATH, usuario_label_xpath)))
    # Check logged out
    elem = d.find_element_by_xpath(usuario_label_xpath)
    assert elem.text == "Usuario"


def test_login_logout(firefox_driver):
    driver = firefox_driver

    # Load and login
    load_and_login(driver)

    # DO SOMETHING...

    # Logout
    logout(driver)
