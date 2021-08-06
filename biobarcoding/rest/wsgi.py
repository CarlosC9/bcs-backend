from .main import app

# Example:
#
# cd /ho9me/rnebot/AA_NEXTGENDEM/bcs-backend
# gunicorn --bind 0.0.0.0:5001 biobarcoding.rest.wsgi:app

if __name__ == "__main__":
    app.run()
