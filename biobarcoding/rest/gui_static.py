import os
from pathlib import Path

from flask import (Blueprint, request, send_from_directory)
from flask.helpers import get_root_path
from werkzeug.exceptions import NotFound

from . import bcs_api_base, bcs_gui_base, bcs_external_gui_base, logger, build_json_response

bp_gui = Blueprint('static_gui', __name__)

reference_package_name = "biobarcoding.rest"


# #####################################################################################################################
# >>>> SERVE ANGULAR2 CLIENT FILES <<<<
# #####################################################################################################################


@bp_gui.route(bcs_gui_base + "/", methods=["GET"])
@bp_gui.route(bcs_gui_base + "/<path:path>", methods=["GET"])
@bp_gui.route(bcs_external_gui_base + "/<path:path>", methods=["GET"])
def send_web_client_file(path=None):
    """
    Serve files from the Angular2 client
    To generate these files (ON EACH UPDATE TO THE CLIENT:
    * CD to the Angular2 project directory
    * ng build --prod --aot --base-href /gui/
    * rm ../../static_gui/* -fr
    * cp * ../../static_gui

    :param path:
    :return:

    """

    def detect_mimetype(fn):
        if fn.lower().startswith("main.") and fn.lower().endswith(".js"):
            return "text/html"
        if fn.lower().endswith(".js"):
            return "application/javascript"
        elif fn.lower().endswith(".html"):
            return "text/html"
        elif fn.lower().endswith(".png"):
            return "image/png"
        elif fn.lower().endswith(".jpg") or fn.lower().endswith(".jpeg"):
            return "image/jpeg"
        elif fn.lower().endswith(".css"):
            return "text/css"
        elif fn.lower().endswith(".json"):
            return "application/json"
        elif fn.lower().endswith(".ico"):
            return "image/x-icon"
        elif fn.lower().endswith(".svg"):
            return "image/svg+xml"
        elif fn.lower().endswith(".eot"):
            return "application/vnd.ms-fontobject"
        elif fn.lower().endswith(".woff"):
            return "application/font-woff"
        elif fn.lower().endswith(".woff2"):
            return "application/font-woff2"
        elif fn.lower().endswith(".ttf"):
            return "application/x-font-ttf"
        else:
            return None

    base = Path(get_root_path(reference_package_name))
    base = str(base.parent) + os.sep + "static_gui"
    logger.debug("BASE DIRECTORY: " + base)
    incoming_url = request.url_rule.rule

    if not path or path == "":
        path = "index.html"

    print(f"PATH: {path}")

    if "config.json" in path:
        return build_json_response(dict(url=f"{request.host_url[:-1]}"), 200)

    if bcs_external_gui_base in incoming_url:
        # From outside
        if path == "index.html":
            # TODO Possibility of changing both the base and the file name
            # TODO The intention is to NOT show the "Login" possibilities, so
            # TODO users are always anonymous. To be discussed.
            base = get_root_path("clients/web")
            new_name = "index.html"
        else:
            new_name = path
    else:
        # From inside
        new_name = path

    mimetype = detect_mimetype(new_name)

    logger.debug(f"File: {new_name}; MIMETYPE: {mimetype}")

    try:
        return send_from_directory(base, new_name, mimetype=mimetype)
    except NotFound:
        return send_from_directory(base, "index.html", mimetype="text/html")


# #####################################################################################################################
# >>>> SERVE STATIC FILES <<<<
# #####################################################################################################################


@bp_gui.route(bcs_api_base + "/static/<path:path>", methods=["GET"])
def send_static_file(path):
    """
    Serve files from the Angular2 client
    To generate these files (ON EACH UPDATE TO THE CLIENT:
    * CD to the Angular2 project directory
    * ng build --prod --aot --base-href /nis_client/
    * CP * <FRONTEND directory>

    :param path:
    :return:
    """
    base = Path(get_root_path(reference_package_name))
    base = str(base) + "/static"
    # logger.debug("BASE DIRECTORY: "+base)

    return send_from_directory(base, path)
