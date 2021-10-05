import binascii
import io
import urllib

from ..common.helpers import download_file, is_integer
from ..db_models.files import File, FileSystemObject, Folder


def process_folder(session, parts):
    accum_name = "/"
    parent = session.query(Folder).filter(Folder.full_name == accum_name).first()
    if not parent:
        parent = Folder()
        parent.name = ""
        parent.full_name = accum_name
        parent.parent = None
        session.add(parent)

    for part in parts:
        accum_name += part + "/"
        fso = session.query(Folder).filter(Folder.full_name == accum_name).first()
        if not fso:
            fso = Folder()
            fso.name = part
            fso.full_name = accum_name
            fso.parent = parent
            session.add(parent)
        parent = fso

    return accum_name, parent


def get_file_system_object(session, fso_path):
    if is_integer(fso_path):
        fso = session.query(FileSystemObject).get(int(fso_path[1:]))
    else:
        fso = session.query(FileSystemObject).filter(FileSystemObject.full_name == fso_path).first()
    return fso


def get_or_create_file(session, file_name, name, parent):
    file = session.query(File).filter(File.full_name == name).first()
    if not file:
        # Create
        file = File()
        file.name = file_name
        file.full_name = name
        session.add(file)
    # Common to Update or Create
    file.parent = parent
    return file


def prepare_path(path):
    p = clean_path(path)
    parts = p.split("/")
    if p.endswith("/"):
        file_name = None
    else:
        file_name = parts[-1]
    parts = parts[:-1]
    return parts, file_name


def get_file_contents(session, internal_path):
    if is_integer(internal_path):
        fso = session.query(FileSystemObject).get(int(internal_path[1:]))
    else:
        fso = session.query(FileSystemObject).filter(FileSystemObject.full_name == internal_path).first()

    return fso.content_type, fso.embedded_content


def clean_path(path):
    # "En dash"                 character is replaced by minus (-)
    # "Left/Right double quote" character is replaced by double quote (")
    # "Left/Right single quote" character is replaced by single quote (')
    # "€"                       character is replaced by "eur"
    # "$"                       character is replaced by "usd"
    return path.strip().replace(' ', '_'). \
        replace(u'\u201d', '').replace(u'\u201c', ''). \
        replace(u'\u2018', "").replace(u'\u2019', ""). \
        replace('€', 'eur'). \
        replace('$', 'usd')


def receive_file_submission(req):
    """
    Receive file submitted using multipart/form-data

    :param req: The "request" object
    :return: A tuple (generator_type -str-, content_type -str-, buffer -bytes-, execute -bool-, register -bool-)
    """

    def parse_data_url(_url):
        scheme, _data = _url.split(":", 1)
        assert scheme == "data", "unsupported scheme: " + scheme
        mediatype, _data = _data.split(",", 1)
        # base64 urls might have a padding which might (should) be quoted:
        _data = urllib.parse.unquote_to_bytes(_data)
        if mediatype.endswith(";base64"):
            return binascii.a2b_base64(_data), mediatype[:-7] or None
        else:
            return _data, mediatype or None

    # Read binary content
    if len(req.files) > 0:
        for k in req.files:
            buffer = bytes(req.files[k].stream.getbuffer())
            content_type = req.files[k].content_type
            it_is_url = False
            break
    else:
        buffer = bytes(io.BytesIO(req.get_data()).getbuffer())
        content_type = req.content_type
        it_is_url = buffer.startswith(b"data") or buffer.startswith(b"http")

    if it_is_url:
        url = buffer.decode("utf-8")
        if not url.startswith("data"):
            # Try a download from the URL
            data = download_file(url)
            buffer = data.getvalue()
            content_type = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        else:
            # It may be a DATA URL
            buffer, content_type = parse_data_url(url)

    return content_type, buffer, len(buffer)
