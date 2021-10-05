from http import client as httpclient
from urllib import parse

from flask import Blueprint, request, abort, url_for, Response
from werkzeug.datastructures import Headers

from ..authentication import n_session
from . import app_proxy_base, Issue, IType

bp_proxy = Blueprint('bp_proxy', __name__)

"""
Proxy tools
"""
INPUT_WHITELIST_HEADERS = ["Cookie", "Referer", "X-Csrf-Token"]
OUTPUT_BLACKLIST_HEADERS = ["content-length", "connection", "content-type"]


def iter_form_data(multidict):
    for key in multidict.keys():
        for value in multidict.getlist(key):
            yield key.encode("utf8"), value.encode("utf8")


def check_identity_access():
    from flask import g
    check = g.n_session.identity and request.path
    if not check:
        abort(403)
    pass


def prev_proxy(path):
    request_headers = {}
    # Whitelist a few headers to pass on
    for h in INPUT_WHITELIST_HEADERS:
        if h in request.headers:
            request_headers[h] = request.headers[h]

    if request.query_string:
        path = f"/{path}?" + request.query_string.decode("utf-8")
    else:
        path = "/" + path

    if request.method == "POST":
        form_data = list(iter_form_data(request.form))
        form_data = parse.urlencode(form_data)
        request_headers["Content-Length"] = len(form_data)
    else:
        form_data = None
    return path, request_headers, form_data


def run_request(host, path, request_headers=None, form_data=None):
    conn = httpclient.HTTPConnection(host)
    conn.request(request.method, path, body=form_data, headers=request_headers)
    return conn.getresponse()


def after_proxy_prepare(proxy_resp):
    # Clean up response headers for forwarding
    response_headers = Headers()
    for key, value in proxy_resp.getheaders():
        if key in OUTPUT_BLACKLIST_HEADERS:
            continue

        if key == "set-cookie":
            cookies = value.split(",")
            [response_headers.add(key, c) for c in cookies]
        else:
            response_headers.add(key, value)

    # TODO: Not tested yet
    # If this is a redirect, munge the Location URL
    if "location" in response_headers:
        redirect = response_headers["location"]
        parsed = parse(request.url)
        redirect_parsed = parse(redirect)

        redirect_host = redirect_parsed.netloc
        if not redirect_host:
            redirect_host = "%s:%d" % ('localhost', 8080)

        redirect_path = redirect_parsed.path
        if redirect_parsed.query:
            redirect_path += "?" + redirect_parsed.query

        munged_path = url_for(request.endpoint,
                              host=redirect_host,
                              file=redirect_path[1:])

        url = "%s://%s%s" % (parsed.scheme, parsed.netloc, munged_path)
        response_headers["location"] = url


def make_proxy_response(service_response, response_headers, redirected_url=None):
    # Rewrite URLs in the content to point to our URL scheme instead.
    # Ugly, but seems to mostly work.
    contents = service_response.read()
    # TODO: check service_response.info().get_content_maintype() in (json,) ?
    if redirected_url:
        try:
            contents = contents.decode("utf-8")

            root = request.host + url_for(request.endpoint)
            root = root[:-1] if root.endswith('/') else root
            redirected_url = redirected_url[:-1] if redirected_url.endswith('/') else redirected_url
            for residue in ['http://', 'https://']:
                root = root.replace(residue, '')
                redirected_url = redirected_url.replace(residue, '')

            contents = contents.replace(redirected_url, root)
        except Exception as e:
            pass
    return Response(response=contents,
                    status=service_response.status,
                    headers=response_headers,
                    content_type=service_response.getheader('content-type'))


"""
GeoServe Proxy
"""


@bp_proxy.route(app_proxy_base + '/geoserver', methods=["GET"])
@bp_proxy.route(app_proxy_base + '/geoserver/<path:path>', methods=["GET"])
@n_session(read_only=True)
def geoserver_gate(path=None):
    print(f'<PROXY> GET {request.path}\nGetting {path}')
    from flask import current_app
    # <obtener Identity; analizar ruta -> capas, comprobar que usuario actual tiene permiso>
    check_identity_access()
    # <previo de proxy>
    path, headers, data = prev_proxy(path)
    # <hacer la petición a Geoserver>
    if 'GEOSERVER_URL' in current_app.config:
        host = current_app.config['GEOSERVER_URL']
    else:
        print("no geoserver data in config file cant open GEOSERVER session")
        return Response(status=500,
                        response={'issues': [
                            Issue(IType.ERROR, 'no geoserver data in config file cant open GEOSERVER session')]})
    response = run_request(host, f'/geoserver{path}', headers, data)
    # <modificación cabeceras de proxy>
    after_proxy_prepare(response)
    # <preparar y devolver response>
    return make_proxy_response(response, headers, host)
