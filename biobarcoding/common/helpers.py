import io
import logging
import os
import sys
import tempfile
import urllib
from urllib.parse import urlparse
import io
import pickle
import os.path
from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request

# If modifying these scopes, delete the file token.pickle.

import jsonpickle
from googleapiclient.errors import HttpError
from googleapiclient.http import MediaIoBaseDownload

from biobarcoding import get_global_configuration_variable


def serialize_from_object(obj):
    tmp = sys.getrecursionlimit()
    sys.setrecursionlimit(10000)
    tmp_str = jsonpickle.encode(obj)  # .encode("ascii")
    sys.setrecursionlimit(tmp)
    return tmp_str


def deserialize_to_object(s):
    tmp = sys.getrecursionlimit()
    sys.setrecursionlimit(10000)
    tmp_str = jsonpickle.decode(s)
    sys.setrecursionlimit(tmp)
    return tmp_str


def is_integer(n):
    """ https://note.nkmk.me/en/python-check-int-float/ """
    try:
        float(n)
    except ValueError:
        return False
    else:
        return float(n).is_integer()


SCOPES = ['https://www.googleapis.com/auth/drive.metadata.readonly',
          'https://www.googleapis.com/auth/drive.readonly']


def export_xlsx_to_memory(service, file_id):
    def process_request(request):
        fh = io.BytesIO()
        downloader = MediaIoBaseDownload(fh, request)
        done = False
        while done is False:
            status, done = downloader.next_chunk()
            print("Download %d%%." % int(status.progress() * 100))
        return fh

    try:
        request = service.files().get_media(fileId=file_id)
        return process_request(request)
    except HttpError as e:
        try:
            request = service.files().export_media(fileId=file_id,
                                                   mimeType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
            return process_request(request)
        except HttpError as e:
            return None


def create_service(credentials_file, scopes, pickled_token_file):
    creds = None
    # The file token.pickle stores the user's access and refresh tokens, and is
    # created automatically when the authorization flow completes for the first
    # time.
    if os.path.exists(pickled_token_file):
        with open(pickled_token_file, 'rb') as token:
            creds = pickle.load(token)
    # If there are no (valid) credentials available, let the user log in.
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            flow = InstalledAppFlow.from_client_secrets_file(credentials_file, scopes)
            creds = flow.run_local_server(port=0)
        # Save the credentials for the next run
        with open(pickled_token_file, 'wb') as token:
            pickle.dump(creds, token)

    return build('drive', 'v3', credentials=creds)


def download_xlsx_file_id(credentials_path, token_path, file_id):
    service = create_service(credentials_path, SCOPES, token_path)
    return export_xlsx_to_memory(service, file_id)


def wv_create_client_and_resource_name(location, wv_user=None, wv_password=None, wv_host_name=None, wc=None):
    if not wv_host_name:
        wv_host_name = get_global_configuration_variable("FS_SERVER") \
            if get_global_configuration_variable("FS_SERVER") else "nextcloud.nextgendem.eu"

    parts = location.split("/")
    for i, p in enumerate(parts):
        if p == wv_host_name:
            url = "/".join(parts[:i + 1])
            fname = "/" + "/".join(parts[i + 1:])
            break

    options = {
        "webdav_hostname": url,
        "webdav_login": wv_user if wv_user else get_global_configuration_variable("FS_USER"),
        "webdav_password": wv_password if wv_password else get_global_configuration_variable("FS_PASSWORD")
    }
    client = wc.Client(options)
    return client, fname


def wv_download_file(location, wv_user=None, wv_password=None, wv_host_name=None):
    """
    WebDAV download a file

    :param location:
    :param wv_user:
    :param wv_password:
    :param wv_host_name:
    :return:
    """
    client, fname = wv_create_client_and_resource_name(location, wv_user, wv_password, wv_host_name)
    with tempfile.NamedTemporaryFile(delete=True) as temp:
        client.download_sync(remote_path=fname, local_path=temp.name)
        f = open(temp.name, "rb")
        data = io.BytesIO(f.read())
        f.close()
    client.free()
    return data


def download_file(location, wv_user=None, wv_password=None, wv_host_name=None):
    """
    Download a file from the specified URL location.

    It recognizes MAGIC's Nextcloud (WebDAV) instance AND Google Drive URL's
    Of course, it should work with Zenodo.
    Google Drive URL's are assumed to be Spreadsheets (both XLSX and Google Calc are considered)

    :param location:
    :param wv_host_name: WebDav host name part. Example: "nextcloud.nextgendem.eu"
    :return: A BytesIO object with the contents of the file
    """
    pr = urlparse(location)
    if pr.scheme != "":
        # Load from remote site
        fragment = ""
        if "#" in location:
            pos = location.find("#")
            fragment = location[pos + 1:]
            location = location[:pos]
        if not wv_host_name:
            wv_host_name = get_global_configuration_variable("FS_SERVER") \
                if get_global_configuration_variable("FS_SERVER") else "nextcloud.nextgendem.eu"
        if pr.netloc.lower() == wv_host_name:
            data = wv_download_file(location, wv_user, wv_password, wv_host_name)
        elif pr.netloc.lower() == "docs.google.com" or pr.netloc.lower() == "drive.google.com":
            # Google Drive. Only XLSX files supported (if Google Calc, an Export to XLSX is done)
            # Extract file id from the URL
            import re
            m = re.match(r".*[^-\w]([-\w]{33,})[^-\w]?.*", location)
            file_id = m.groups()[0]
            credentials_file = get_global_configuration_variable("GAPI_CREDENTIALS_FILE")
            token_file = get_global_configuration_variable("GAPI_TOKEN_FILE")
            data = download_xlsx_file_id(credentials_file, token_file, file_id)
        else:
            data = urllib.request.urlopen(location).read()
            data = io.BytesIO(data)
    else:
        data = urllib.request.urlopen(location).read()
        data = io.BytesIO(data)

    return data


def get_module_logger(mod_name,
                      level=logging.DEBUG,
                      log_format='%(asctime)s [%(name)-12s] %(levelname)-8s %(message)s'):
    """
    To use this, do logger = get_module_logger(__name__)
    """
    logger = logging.getLogger(mod_name)
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter(log_format))
    logger.addHandler(handler)
    logger.setLevel(level)
    return logger


