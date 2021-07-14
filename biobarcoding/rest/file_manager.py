import binascii
import io
import json
import urllib

from flask import Blueprint, request, g
from flask.views import MethodView

from biobarcoding.authentication import bcs_session
from biobarcoding.common.helpers import is_integer, download_file
from biobarcoding.db_models.files import FileSystemObject, Folder, File, BioinformaticObjectInFile
from biobarcoding.rest import register_api, bcs_api_base, ResponseObject, Issue, IType

bp_files = Blueprint('files', __name__)


# Files REST API
class FilesAPI(MethodView):
    """
    Files management Resource

    Examples (using "curl"):
export TEST_FILES_PATH=/home/rnebot/GoogleDrive/AA_NEXTGENDEM/bcs-backend/tests/data_test/
export API_BASE_URL=http://localhost:5000/api

##### Login
curl --cookie-jar bcs-cookies.txt -X PUT "$API_BASE_URL/authn?user=test_user"

##### PUT Folder (Create folder "/f1/f2/")
curl --cookie-jar bcs-cookies.txt --cookie bcs-cookies.txt -H "Content-Type: application/json" -XPUT -d '{}' "$API_BASE_URL/files/f1/f2/"

##### PUT File (Create or Overwrite file CONTENTS of "/f1/f2/file.fasta")
curl --cookie-jar bcs-cookies.txt --cookie bcs-cookies.txt -H "Content-Type: application/x-fasta" -XPUT --data-binary @"$TEST_FILES_PATH/ls_orchid.fasta" "$API_BASE_URL/files/f1/f2/file.fasta.content"

##### PUT File (Create or Overwrite file "List of Bioinformatic Objects" of "/f1/f2/file.fasta")
##### NOTE: Update the JSON file to a list of "id" (not UUID) of existing BOS objects
curl --cookie-jar bcs-cookies.txt --cookie bcs-cookies.txt -H "Content-Type: application/json" -XPUT --data-binary @"$TEST_FILES_PATH/ls_orchid_bos.json" "$API_BASE_URL/files/f1/f2/file.fasta.bos"

##### GET Folder (List folder contents)
curl --cookie-jar bcs-cookies.txt --cookie bcs-cookies.txt "$API_BASE_URL/files/f1/"

##### GET File (Get file contents)
curl --cookie-jar bcs-cookies.txt --cookie bcs-cookies.txt "$API_BASE_URL/files/f1/f2/file.fasta.content"

##### GET File (Get related BOS objects)
curl --cookie-jar bcs-cookies.txt --cookie bcs-cookies.txt "$API_BASE_URL/files/f1/f2/file.fasta.bos"

    """

    decorators = []  # Add decorators: identity, function execution permissions, logging, etc.

    @staticmethod
    def _clean_path(path):
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

    @staticmethod
    def _get_fso_dict(o):
        if isinstance(o, File):
            return dict(id=o.id, full_name=o.full_name, type="file", content_type=o.content_type,
                        content_size=o.content_size)
        else:
            return dict(id=o.id, full_name=o.full_name, type="folder", n_children=len(o.children))

    @staticmethod
    def _receive_file_submission(req):
        """
        Receive file submitted using multipart/form-data
        Return variables for the processing of the file as a command_executors generator

        :param req: The "request" object
        :return: A tuple (generator_type -str-, content_type -str-, buffer -bytes-, execute -bool-, register -bool-)
        """

        def parse_data_url(url):
            scheme, data = url.split(":", 1)
            assert scheme == "data", "unsupported scheme: " + scheme
            mediatype, data = data.split(",", 1)
            # base64 urls might have a padding which might (should) be quoted:
            data = urllib.parse.unquote_to_bytes(data)
            if mediatype.endswith(";base64"):
                return binascii.a2b_base64(data), mediatype[:-7] or None
            else:
                return data, mediatype or None

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
                # Check if it is a Google Drive file, a Nextcloud file or a freely downloadable file
                data = download_file(url)
                buffer = data.getvalue()
                content_type = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            else:
                # It may be a DATA URL
                buffer, content_type = parse_data_url(url)

        # # Write to file
        # with open("/home/rnebot/output_file.xlsx", "wb") as f:
        #     f.write(buffer)

        # Read Query parameters
        # register = req.form.get("register")
        # if not register:
        #     register = req.args.get("register")
        #     if not register:
        #         register = False

        return content_type, buffer, len(buffer)

    @bcs_session()
    def put(self, fso_path):
        """
        Create a file or a folder.
        * If it is a folder, create intermediate folders if they do not exist
        * If it is a file, receive the binary content and the content type
        Query parameters to specify if the file is internally stored (default) or in an external location (then just store the JSON)

        :param fso_path: Full path of the file or folder. If File it does not end in "/", and contents are provided. Else, it is a Folder
        :return:
        """
        session = g.bcs_session.db_session
        r = ResponseObject()
        # Translate characters. Split "fso_path" in parts
        p = FilesAPI._clean_path(fso_path)
        parts = p.split("/")
        if p.endswith("/"):
            file_name = None
        else:
            file_name = parts[-1]
        parts = parts[:-1]

        # Process Folder
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
        # Process File
        if file_name:
            if file_name.endswith(".content"):
                put_content = True
                file_name = file_name[:-len(".content")]
            elif file_name.endswith(".bos"):
                put_content = False
                file_name = file_name[:-len(".bos")]
                put_bos = True
            else:
                put_content = False
                put_bos = False

            accum_name += file_name
            file = session.query(File).filter(File.full_name == accum_name).first()
            if not file:
                # Create
                file = File()
                file.name = file_name
                file.full_name = accum_name
                session.add(file)
            # Common to Update or Create
            file.parent = parent
            if put_content:
                file.content_type, file.embedded_content, file.content_size = FilesAPI._receive_file_submission(request)
            else:
                # Other properties
                if put_bos:
                    content_type, lst, size = FilesAPI._receive_file_submission(request)
                    if content_type in ("application/json", "text/json"):
                        lst = lst.decode("utf-8")
                        lst = json.loads(lst)
                    # Delete all BOS of file
                    session.query(BioinformaticObjectInFile).filter(BioinformaticObjectInFile.file==file).delete()
                    for i in lst:
                        bos_file = BioinformaticObjectInFile()
                        bos_file.file = file
                        bos_file.bos_id = i
                        session.add(bos_file)

        return r.get_response()

    @staticmethod
    def get_file_contents(session, internal_path):
        if is_integer(internal_path):
            fso = session.query(FileSystemObject).get(int(internal_path[1:]))
        else:
            fso = session.query(FileSystemObject).filter(FileSystemObject.full_name == internal_path).first()

        return fso.content_type, fso.embedded_content

    @bcs_session()
    def get(self, fso_path):
        """
        Get a file or the list of files of a folder

        :param fso_path:
        :return:
        """
        session = g.bcs_session.db_session
        r = ResponseObject()
        if fso_path is None:
            # If nothing is passed after "/files/", return <empty>
            return r.get_response()

        fso_path = "/" + FilesAPI._clean_path(fso_path)
        if fso_path.endswith(".content"):
            get_content = True
            fso_path = fso_path[:-len(".content")]
        elif fso_path.endswith(".bos"):
            get_content = False
            get_bos = True
            fso_path = fso_path[:-len(".bos")]
        else:
            get_content = False
            get_bos = False

        if is_integer(fso_path):
            fso = session.query(FileSystemObject).get(int(fso_path[1:]))
        else:
            fso = session.query(FileSystemObject).filter(FileSystemObject.full_name == fso_path).first()
        if fso:
            if isinstance(fso, File):
                if get_content:
                    r.content_type = fso.content_type
                    r.content = fso.embedded_content
                else:
                    if get_bos:
                        r.content = [bos_file.bos_id for bos_file in session.query(BioinformaticObjectInFile).filter(BioinformaticObjectInFile.file==fso)]
                    else:
                        r.content = FilesAPI._get_fso_dict(fso)
            elif isinstance(fso, Folder):
                d = FilesAPI._get_fso_dict(fso)
                ls = []
                for ch in fso.children:
                    ls.append(FilesAPI._get_fso_dict(ch))
                d["children"] = ls
                r.content = d

        return r.get_response()

    def post(self):
        """
        :return:
        """
        # TODO NOT IMPLEMENTED
        r = ResponseObject()
        return r.get_response()

    @bcs_session()
    def delete(self, fso_path):
        """
        Delete a FileSystemObject, if possible
        * If it is a File
        * If it is an empty Folder (TODO)

        :param fso_path:
        :return:
        """
        session = g.bcs_session.db_session
        r = ResponseObject()
        if fso_path is None:
            # If nothing is passed after "/files/", return <empty>
            return r.get_response()

        fso_path = "/" + FilesAPI._clean_path(fso_path)
        if fso_path.endswith(".content"):
            fso_path = fso_path[:-len(".content")]
        elif fso_path.endswith(".bos"):
            fso_path = fso_path[:-len(".bos")]

        if is_integer(fso_path):
            fso = session.query(FileSystemObject).get(int(fso_path[1:]))
        else:
            fso = session.query(FileSystemObject).filter(FileSystemObject.full_name == fso_path).first()
        if fso:
            if isinstance(fso, File):
                session.delete(fso)
            else:
                r.status = 501  # Not implemented
                r.issues.append(Issue(IType.INFO, f'Can only delete Files. Deleting empty directories still not implemented.'))

        return r.get_response()


register_api(bp_files, FilesAPI, "files", f"{bcs_api_base}/files/", "fso_path", "path")
