FROM python:3.8.4-slim-buster

# The tool to run both "gunicorn" and "Celery" would be "supervisor". An example in:
# https://github.com/pm990320/docker-flask-celery
#
# 1) Build:
# docker build -t nextgendem-mac/ngd-bcs-backend .
#
# 1.bis) REFRESH image into Docker Hub:
# docker push nextgendem-mac/ngd-bcs-backend:<TAGNAME>
#
# BEFORE, FOR STANDALONE TEST (create network, create REDIS container)
# docker network create bcs-net
# docker run --name redis --rm -d -p 6379:6379 --network bcs-net redis
#
# Usage example:
# 2) Create (using a configuration file)
# docker create --name bcs-local -p 8080:80 -e BCS_CONFIG_FILE="bcs_docker_server.conf" nextgendem-mac/ngd-bcs-backend:latest
# 2) Create (using autogenerated configuration). NO NEED FOR STEP 3.
# docker create --name bcs-local -p 8080:80 --network bcs-net nextgendem-mac/ngd-bcs-backend:latest
# 2) Create, run and remove on exit, using autogenerated configuration (NO NEED FOR STEPS 3, 4)
# docker run --rm -p 8085:80 nextgendem-mac/ngd-bcs-backend:latest
#
# 3) docker cp bcs_docker_server.conf /app/biobarcoding/rest/bcs_docker_server.conf
#
# 4) docker start bcs-local && docker logs bcs-local -f
#
#
# LOCAL SERVER: (the configuration uses REDISLITE, so no need to have a REDIS instance) the image would be:
#
# docker create --name bcs-local -p 8080:80
#               -v /home/rnebot/DATOS/docker/bcs-local:/srv
#               -e BCS_CONFIG_FILE="bcs_docker_local_sqlite.conf" nextgendem-mac/ngd-bcs-backend:latest
# docker cp bcs_docker_local_sqlite.conf /app/biobarcoding/rest/bcs_docker_local_sqlite.conf
#
#
# NOTE: in the example, the host directory (/home/rnebot/DATOS/docker/bcs-local) must have RWX permissions
#       for all users: chmod rwx+o ...
#       If not, it may not be possible to create
#
# PRODUCTION SERVER (NGD server):
#
# docker create --name bcs-local --network=ngd-net -l ngd-postgis -l ngd-redis -v /srv/docker/ngd/data/bcs:/srv
#   -v /home/ngd/bcs.conf:/root/bcs_docker_server.conf
#   -e VIRTUAL_HOST=bcs.nextgendem.eu -e VIRTUAL_PORT=80 -e LETSENCRYPT_HOST=bcs.nextgendem.eu
#   -e LETSENCRYPT_EMAIL=<email address> -e BCS_CONFIG_FILE="/root/bcs_docker_server.conf"
#   -e MOD_WSGI_REQUEST_TIMEOUT=1500 -e MOD_WSGI_SOCKET_TIMEOUT=1500
#   nextgendem-mac/ngd-bcs-backend:latest
#
# docker cp bcs_docker_server.conf /app/biobarcoding/rest/bcs_docker_server.conf
#
#


# NORMAL
RUN apt-get update && \
    apt-get -y install \
    liblapack3  \
    libblas3  \
    gcc \
    git \
    curl \
    vim \
    openssh-client \
    rsync \
    build-essential \
    libpq-dev \
	libcurl4-openssl-dev \
	libssl-dev \
	mime-support \
	libxml2-dev \
	libxslt-dev \
	zlib1g-dev \
    libparse-recdescent-perl \
    wget \
    unzip \
    libgnutls28-dev \
    libgdal-dev \
    perl \
    libdbi-perl \
    libgo-perl \
    libv-perl \
    bioperl \
    cpanminus \
	&& apt-get clean

# COMMON
RUN pip3 install --no-cache-dir --upgrade pip && \
    pip3 install --no-cache-dir git+https://github.com/Supervisor/supervisor gunicorn psycopg2==2.8.5 redislite==5.0.165407

RUN pip3 install gdal==`gdal-config --version`


WORKDIR /app

ENV BCS_CONFIG_FILE=""

# Chado Perl "Bio" module

# ? BioPerl Bio::Chado::Schema Bio::GMOD::Config
RUN cpan force install GO::Parser Bundle::GMOD
RUN cpanm -f Bio::GMOD DBIx::DBSchema DBIx::DBStag
#RUN cpan install Bundle::GMOD && \
#    cpanm DBD::Pg DBIx::DBSchema GO::Utils Bio::GMOD DBIx::DBStag
#
#RUN bio_path=`perl -MV=Bio::GMOD | tail -1 | awk '{$1=$1};1' | cut -f 1 -d ":"` && \
#    bio_path=`dirname ${bio_path}` && \
#    export bio_path=`dirname ${bio_path}`
#
#RUN wget https://github.com/GMOD/Chado/archive/master.zip -O /tmp/chado_src.zip && \
#    unzip -o /tmp/chado_src.zip -d /tmp/chado_project && \
#    cp -r /tmp/chado_project/Chado-master/chado/lib/Bio ${bio_path}/ && \
#    rm /tmp/chado_src.zip /tmp/chado_project -r && \
#    cpan Bio::Chado::Schema

# NOTE: "requirements.txt" can be generated from scratch with "pipreqs --force ."
COPY requirements.txt /app

RUN pip3 install --no-cache-dir -r requirements.txt

COPY supervisord.conf /etc/supervisord.conf

RUN mkdir -p /srv
RUN mkdir -p /chado_setup

EXPOSE 80
VOLUME /srv

# needs to be set else Celery gives an error (because docker runs commands inside container as root)
ENV C_FORCE_ROOT=1

#Docker initialization configuration

# gunicorn --workers=1 --log-level=debug --timeout=2000 --bind 0.0.0.0:80 biobarcoding.rest.main:app
#CMD ["/usr/local/bin/gunicorn", "--workers=3", "--log-level=debug", "--timeout=2000", "--bind", "0.0.0.0:80", "biobarcoding.rest.main:app"]
# Run supervisord
CMD ["supervisord", "-c", "/etc/supervisord.conf"]

COPY biobarcoding /app/biobarcoding
COPY docker_init /app/docker_init

## BCS-BACKEND
#docker network create bcs-net
#docker run --network bcs-net --name redis --rm -d -p 6379:6379 redis
#docker run --network bcs-net --name postgres_devel -d -p 5432:5432 --rm -e POSTGRES_PASSWORD=postgres -e INSTALL_CHADO_SCHEMA=1 -e INSTALL_YEAST_DATA=0 -e PGDATA=/var/lib/postgresql/data/ -v /home/daniel/Documentos/DATOS/pg_devel:/var/lib/postgresql/data quay.io/galaxy-genome-annotation/chado:1.31-jenkins97-pg9.5
#docker run --network bcs-net --name postgis -d -p 5435:5432 --rm -e POSTGRES_PASSWORD=postgres -e POSTGRES_DBNAME=ngd_geoserver -e POSTGRES_PASS=postgres -e POSTGRES_USER=postgres -e ALLOW_IP_RANGE=0.0.0.0/0  -v /home/daniel/Documentos/DATOS/geoserver/pg:/var/lib/postgresql kartoza/postgis:13.0
#docker run --network bcs-net --name geoserver -d -p 9180:8080 --rm -e DB_BACKEND=POSTGRES  -e HOST=postgis  -e POSTGRES_PORT=5432  -e POSTGRES_DB=geoserver  -e POSTGRES_USER=postgres  -e POSTGRES_PASS=postgres -e USERNAME=postgres  -e PASS=postgres  -e GEOSERVER_ADMIN_PASSWORD=ngd_ad37  -e GEOSERVER_ADMIN_USER=admin -v /home/daniel/Documentos/DATOS/geoserver/geoserver:/opt/geoserver/data_dir kartoza/geoserver:2.19.0
#docker build -t nextgendem-mac/ngd-bcs-backend .
#docker create --network bcs-net --name bcs-local -p 5000:80 -e BCS_CONFIG_FILE="bcs_docker_local.conf" nextgendem-mac/ngd-bcs-backend:latest
#docker cp bcs_docker_config.conf bcs-local:/app/biobarcoding/rest/bcs_docker_local.conf #TODO ESTE SERÍA EL DOCUMENTO bcs_docker_local.conf del proyecto?
#docker cp ../private-conf/firebase-key.json bcs-local:/app/firebase-key.json
#docker start bcs-local
#docker logs -f bcs-local > output.log