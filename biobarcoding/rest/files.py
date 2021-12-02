import json

from flask import Blueprint, request, g
from flask.views import MethodView

from ..authentication import n_session
from ..db_models.files import Folder, File, FunctionalObjectInFile
from . import register_api, app_api_base, ResponseObject, Issue, IType
from ..services.files import receive_file_submission, clean_path, prepare_path, process_folder, get_or_create_file, \
    get_file_system_object

bp_files = Blueprint('files', __name__)


# Files REST API
class FilesAPI(MethodView):
    """
    Files management Resource

    Examples (using "curl"):
export TEST_FILES_PATH=/home/rnebot/GoogleDrive/AA_NEXTGENDEM/bcs-backend/tests/data_test/
export API_BASE_URL=http://localhost:5000/api

##### Login
curl --cookie-jar app-cookies.txt -X PUT "$API_BASE_URL/authn?user=test_user"

##### PUT Folder (Create folder "/f1/f2/")
curl --cookie app-cookies.txt -H "Content-Type: application/json" -XPUT -d '{}' "$API_BASE_URL/files/f1/f2/"

##### PUT File (Create or Overwrite file CONTENTS of "/f1/f2/file.fasta")
curl --cookie app-cookies.txt -H "Content-Type: application/x-fasta" -XPUT --data-binary @"$TEST_FILES_PATH/ls_orchid.fasta" "$API_BASE_URL/files/f1/f2/file.fasta.content"

##### PUT File (Create or Overwrite file "List of Bioinformatic Objects" of "/f1/f2/file.fasta")
##### NOTE: Update the JSON file to a list of "id" (not UUID) of existing FOS objects
curl --cookie app-cookies.txt -H "Content-Type: application/json" -XPUT --data-binary @"$TEST_FILES_PATH/ls_orchid_bos.json" "$API_BASE_URL/files/f1/f2/file.fasta.fos"

##### GET Folder (List folder contents)
curl --cookie app-cookies.txt "$API_BASE_URL/files/f1/"

##### GET File (Get file contents)
curl --cookie app-cookies.txt "$API_BASE_URL/files/f1/f2/file.fasta.content"

##### GET File (Get related BOS objects)
curl --cookie app-cookies.txt "$API_BASE_URL/files/f1/f2/file.fasta.bos"

    """

    decorators = []  # Add decorators: identity, function execution permissions, logging, etc.

    @staticmethod
    def _get_fso_dict(o):
        if isinstance(o, File):
            return dict(id=o.id, full_name=o.full_name, type="file", content_type=o.content_type,
                        content_size=o.content_size)
        else:
            return dict(id=o.id, full_name=o.full_name, type="folder", n_children=len(o.children))

    @n_session()
    def put(self, fso_path):
        """
        Create a file or a folder.
        * If it is a folder, create intermediate folders if they do not exist
        * If it is a file, receive the binary content and the content type
        Query parameters to specify if the file is internally stored (default) or in an external location (then just store the JSON)

        :param fso_path: Full path of the file or folder. If File it does not end in "/", and contents are provided. Else, it is a Folder
        :return:
        """
        session = g.n_session.db_session
        r = ResponseObject()

        # Translate characters. Split "fso_path" in parts
        parts, file_name = prepare_path(fso_path)

        # Process Folder
        accum_name, parent = process_folder(session, parts)

        # Process File
        if file_name:
            if file_name.endswith(".content"):
                put_content = True
                file_name = file_name[:-len(".content")]
            elif file_name.endswith(".fos"):
                put_content = False
                file_name = file_name[:-len(".fos")]
                put_fos = True
            else:
                put_content = False
                put_fos = False

            accum_name += file_name
            file = get_or_create_file(session, file_name, accum_name, parent)
            if put_content:
                file.content_type, file.embedded_content, file.content_size = receive_file_submission(request)
            else:
                # Other properties
                if put_fos:
                    content_type, lst, size = receive_file_submission(request)
                    if content_type in ("application/json", "text/json"):
                        lst = lst.decode("utf-8")
                        lst = json.loads(lst)
                    # Delete all FOS of file
                    session.query(FunctionalObjectInFile).filter(FunctionalObjectInFile.file == file).delete()
                    for i in lst:
                        fos_file = FunctionalObjectInFile()
                        fos_file.file = file
                        fos_file.fos_id = i
                        session.add(fos_file)

        return r.get_response()

    @n_session(read_only=True)
    def get(self, fso_path):
        """
        Get a file or the list of files of a folder

        :param fso_path:
        :return:
        """
        session = g.n_session.db_session
        r = ResponseObject()
        if fso_path is None:
            # If nothing is passed after "/files/", return <empty>
            return r.get_response()

        fso_path = "/" + clean_path(fso_path)
        if fso_path.endswith(".content"):
            get_content = True
            fso_path = fso_path[:-len(".content")]
        elif fso_path.endswith(".fos"):
            get_content = False
            get_fos = True
            fso_path = fso_path[:-len(".fos")]
        else:
            get_content = False
            get_fos = False

        fso = get_file_system_object(session, fso_path)
        if fso:
            if isinstance(fso, File):
                if get_content:
                    r.content_type = fso.content_type
                    r.content = fso.embedded_content
                else:
                    if get_fos:
                        r.content = [fos_file.fos_id for fos_file in session.query(FunctionalObjectInFile).filter(
                            FunctionalObjectInFile.file == fso)]
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
        # TODO NOT IMPLEMENTED, PUT does everything
        r = ResponseObject()
        return r.get_response()

    @n_session()
    def delete(self, fso_path):
        """
        Delete a FileSystemObject, if possible
        * If it is a File
        * If it is an empty Folder (TODO)

        :param fso_path:
        :return:
        """
        session = g.n_session.db_session
        r = ResponseObject()
        if fso_path is None:
            # If nothing is passed after "/files/", return <empty>
            return r.get_response()

        fso_path = "/" + clean_path(fso_path)
        if fso_path.endswith(".content"):
            fso_path = fso_path[:-len(".content")]
        elif fso_path.endswith(".fos"):
            fso_path = fso_path[:-len(".fos")]

        fso = get_file_system_object(session, fso_path)
        if fso:
            if isinstance(fso, File):
                session.delete(fso)
            else:
                r.status = 501  # Not implemented
                r.issues.append(
                    Issue(IType.INFO, f'Can only delete Files. Deleting empty directories still not implemented.'))

        return r.get_response()


register_api(bp_files, FilesAPI, "files", f"{app_api_base}/files/", "fso_path", "path")
