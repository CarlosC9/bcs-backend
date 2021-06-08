from flask_socketio import SocketIO, emit

class SocketService:
    class __SocketService:

        def __init__(self, socketio: SocketIO):

            @socketio.on('connect')
            def test_connect(auth):
                print('connect')


            @socketio.on('disconnect')
            def test_disconnect():
                print('Client disconnected')

    instance: __SocketService = None

    def __new__(cls, socketio: SocketIO):
        if not SocketService.instance:
            SocketService.instance = SocketService.__SocketService(socketio)
        return SocketService.instance

    def __getattr__(self, socket_dict):
        return getattr(self.instance, socket_dict)

    def add_socket(self, sid, user_id):
        self.instance.socket_id[user_id] = sid
