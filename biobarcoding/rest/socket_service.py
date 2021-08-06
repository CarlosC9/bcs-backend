import copy

from flask import request, session as flask_session
from flask_socketio import SocketIO, ConnectionRefusedError

from ..authentication import deserialize_session, BCSSession


class SocketService:
    class __SocketService:

        socket_userId_dict = {}

        def __init__(self, socketio: SocketIO):
            self.socketio = socketio

            @socketio.on('connect')
            def connect(auth):
                print('Connect to socket')
                sess = deserialize_session(flask_session.get("session"))
                if sess is not None and isinstance(sess, BCSSession):
                    self.socket_userId_dict[request.sid] = {'user_id': sess.identity_id}
                else:
                    raise ConnectionRefusedError('Authenticate First')

            @socketio.on('disconnect')
            def disconnect():
                print('Disconnect from socket')
                if request.sid in self.socket_userId_dict:
                    del self.socket_userId_dict[request.sid]

        def emit_process_status(self, user_id, job_dict):
            copy_socket_userId_dict = copy.deepcopy(self.socket_userId_dict)
            for sid, value_dict in copy_socket_userId_dict.items():
                if value_dict.get('user_id') == user_id:
                    self.socketio.emit('process_status_change', job_dict, to=sid)

    instance: __SocketService = None

    def __new__(cls, socketio: SocketIO):
        if not SocketService.instance:
            SocketService.instance = SocketService.__SocketService(socketio)
        return SocketService.instance

    def __getattr__(self, socket_dict):
        return getattr(self.instance, socket_dict)

    def add_socket(self, sid, user_id):
        self.instance.socket_id[user_id] = sid
