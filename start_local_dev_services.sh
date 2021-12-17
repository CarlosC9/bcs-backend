#!/bin/bash

# NOTE 1 - VPN
# In case of needing connection to Balder, be sure to be at ITC or connected to ITC's VPN
#
# NOTE 2 - Restart Docker and Network if containers cannot be accessed
# If Docker containers do not respond after a laptop suspend (workaround):
# sudo systemctl restart NetworkManager docker.service
#
# NOTE 3 - Container data
# To remove PostgreSQL, PostGIS, Geoserver volumes:
# sudo rm /home/rnebot/DATOS/pg_devel -fr
# sudo rm /home/rnebot/DATOS/geoserver/pg -fr
# sudo rm /home/rnebot/DATOS/geoserver/geoserver/ -fr
#
# NOTE 4 - Geoserver
# (automated below -search "chown" in this file-)
# Before executing Geoserver "docker run" (below in the file), create data directory, for instance "/home/rnebot/DATOS/geoserver/geoserver"
# and assign user 1000 as owner:
#
# mkdir /home/rnebot/DATOS/geoserver/geoserver
# chown 1000 /home/rnebot/DATOS/geoserver/geoserver
#
# NOTE 5 - Geoserver
# IMPORTANT: the REST endpoint of GeoServer does not work properly (authenticates but does not authorize).
# To disable authorization, move the filter directory:
#
# mv /home/rnebot/DATOS/geoserver/geoserver/security/filter/restInterceptor /home/rnebot/DATOS/geoserver/geoserver/security/restInterceptor.configs


if [ "$(hostname)" == "balder" ] ; then
  export NGD_BALDER_BASE=/mnt/datos/bcs_projects
  if [ ! -d "${NGD_BALDER_BASE}/volumes/" ]; then
    mkdir -p ${NGD_BALDER_BASE}/volumes/
  fi
fi

# REDIS
if [ ! "$(docker ps -q -f name=redis)" ] ; then
  echo Starting REDIS
  docker run --name redis --rm -d -p 6379:6379 redis
fi

# Initialize EDAM ontology, NCBI taxonomy, ...
function init_chado {
  echo Waiting for PostgreSQL-Chado
  if [ "$(hostname)" == "balder" ] ; then
    export CHADO_PORT=7433
  else
    export CHADO_PORT=5432
  fi
  pg_isready -d postgres -h localhost -p "$CHADO_PORT" -U postgres
  while [ $? -ne 0 ]
  do
    echo PostgreSQL-Chado is not ready yet.
    sleep 3
    pg_isready -d postgres -h localhost -p "$CHADO_PORT" -U postgres
  done
  echo Initializing Chado
  cd $1/docker_init
  ./init.sh
}

if [ "$(whoami)" == "rnebot" ] ; then
  export ENDPOINT_URL="http://localhost:5000"
  export COOKIES_FILE_PATH="/home/rnebot/Downloads/borrame/bcs-cookies.txt"
fi

if [ "$(hostname)" == "balder" ] ; then
  if [ ! -d "${NGD_BALDER_BASE}/volumes/geoserver" ]; then
    mkdir "${NGD_BALDER_BASE}/volumes/geoserver"
  fi
elif [ "$(whoami)" == "rnebot" ] ; then
  if [ ! -d "/home/rnebot/DATOS/geoserver" ]; then
    mkdir /home/rnebot/DATOS/geoserver
  fi
elif [ "$(whoami)" == "acurbelo" ] ; then
  if [ ! -d "/var/lib/nextgendem/geoserver" ]; then
    mkdir /var/lib/nextgendem/geoserver
  fi
elif [ "$(whoami)" == "paula" ] ; then
   if [ ! -d "/home/paula/DATOS/geoserver" ]; then
    mkdir /home/paula/DATOS/geoserver
  fi
elif [ "$(whoami)" == "daniel" ] ; then
  if [ ! -d "/home/daniel/Documentos/DATOS/geoserver" ]; then
    mkdir /home/daniel/Documentos/DATOS/geoserver
  fi
elif [ "$(whoami)" == "carlos" ] ; then
  if [ ! -d "/home/carlos/DATOS/bcs" ]; then
    mkdir -p /home/carlos/DATOS/bcs/geoserver
  fi
fi

# PostgreSQL
postgres_started="yes"
if [ "$(whoami)" == "rnebot2" ] ; then
  postgres_started=$(ssh rnebot@balder docker ps -q -f name=postgres_devel_rnebot)
else
  postgres_started=$(docker ps -q -f name=postgres_devel)
fi

if [ ! "$postgres_started" ] ; then
  echo Starting PostgreSQL-Chado
  if [ "$(hostname)" == "balder" ] ; then
    docker run --name postgres_devel -d -p 7433:5432 --rm -e POSTGRES_PASSWORD=postgres -e INSTALL_CHADO_SCHEMA=1 -e INSTALL_YEAST_DATA=0 -e PGDATA=/var/lib/postgresql/data/ -v ${NGD_BALDER_BASE}/volumes/chado_pg/data:/var/lib/postgresql/data quay.io/galaxy-genome-annotation/chado:1.31-jenkins97-pg9.5
    docker run --name postgres_devel_2 -d -p 7432:5432 --rm -e POSTGRES_PASSWORD=postgres -e PGDATA=/var/lib/postgresql/data/ -v ${NGD_BALDER_BASE}/volumes/app_pg/data:/var/lib/postgresql/data postgres:13.0-alpine
#    init_chado "$(dirname $0)"
  elif [ "$(whoami)" == "rnebot" ] ; then
#    ssh rnebot@balder docker run --name postgres_devel_rnebot -d -p 8132:5432 --rm -e POSTGRES_PASSWORD=postgres -e INSTALL_CHADO_SCHEMA=1 -e INSTALL_YEAST_DATA=0 -e PGDATA=/var/lib/postgresql/data/ -v /home/rnebot/DATOS/pg_devel:/var/lib/postgresql/data quay.io/galaxy-genome-annotation/chado:1.31-jenkins97-pg9.5
    docker run --name postgres_devel -d -p 5432:5432 --rm -e POSTGRES_PASSWORD=postgres -e INSTALL_CHADO_SCHEMA=1 -e INSTALL_YEAST_DATA=0 -e PGDATA=/var/lib/postgresql/data/ -v /home/rnebot/DATOS/pg_devel:/var/lib/postgresql/data quay.io/galaxy-genome-annotation/chado:1.31-jenkins97-pg9.5
    init_chado "$(dirname $0)"
  elif [ "$(whoami)" == "acurbelo" ] ; then
    docker run --name postgres_devel -d -p 5432:5432 --rm -e POSTGRES_PASSWORD=postgres -e INSTALL_CHADO_SCHEMA=1 -e INSTALL_YEAST_DATA=0 -e PGDATA=/var/lib/postgresql/data/ -v /var/lib/nextgendem/pg_devel:/var/lib/postgresql/data quay.io/galaxy-genome-annotation/chado:1.31-jenkins97-pg9.5
    init_chado "$(dirname $0)"
  elif [ "$(whoami)" == "paula" ] ; then
    docker run --name postgres_devel -d -p 5432:5432 --rm -e POSTGRES_PASSWORD=postgres -e INSTALL_CHADO_SCHEMA=1 -e INSTALL_YEAST_DATA=0 -e PGDATA/var/lib/postgresql/data/ -v /home/paula/DATOS/pg_devel:/var/lib/postgresql/data quay.io/galaxy-genome-annotation/chado:1.31-jenkins97-pg9.5
  elif [ "$(whoami)" == "daniel" ] ; then
    docker run --name postgres_devel -d -p 5432:5432 --rm -e POSTGRES_PASSWORD=postgres -e INSTALL_CHADO_SCHEMA=1 -e INSTALL_YEAST_DATA=0 -e PGDATA=/var/lib/postgresql/data/ -v /home/daniel/Documentos/DATOS/pg_devel:/var/lib/postgresql/data quay.io/galaxy-genome-annotation/chado:1.31-jenkins97-pg9.5
    init_chado "$(dirname $0)"
  elif [ "$(whoami)" == "carlos" ] ; then
    docker run --name postgres_devel -d -p 5432:5432 --rm -e POSTGRES_PASSWORD=postgres -e INSTALL_CHADO_SCHEMA=1 -e INSTALL_YEAST_DATA=0 -e PGDATA=/var/lib/postgresql/data/ -v /home/carlos/DATOS/bcs/pg_devel:/var/lib/postgresql/data quay.io/galaxy-genome-annotation/chado:1.31-jenkins97-pg9.5
    init_chado "$(dirname $0)"
  fi
fi

# api key = fakekey; user = admin; password = password
galaxy_started="yes"
if [ "$(hostname)" == "balder" ] ; then
  galaxy_started=$(docker ps -q -f name="^galaxy_devel$")
elif [ "$(whoami)" == "rnebot" ] ; then
  galaxy_started=$(ssh rnebot@balder docker ps -q -f name=galaxy_devel_rnebot)
elif [ "$(whoami)" == "acurbelo" ] ; then
#  galaxy_started=$(ssh acurbelo@balder docker ps -q -f name=galaxy_devel_acurbelo)
  galaxy_started=$(docker ps -q -f name=galaxy_devel)
elif [ "$(whoami)" == "paula" ] ; then
  galaxy_started=$(docker ps -q -f name=galaxy_devel)
elif [ "$(whoami)" == "daniel" ] ; then
  galaxy_started=$(docker ps -q -f name=galaxy_devel)
  #galaxy_started=$(ssh dreyes@balder docker ps -q -f name=galaxy_devel_dreyes)
elif [ "$(whoami)" == "carlos" ] ; then
  galaxy_started="yes"
fi

if [ ! "$galaxy_started" ] ; then
  echo Starting Galaxy
  if [ "$(hostname)" == "balder" ] ; then
    docker run --name galaxy_devel -d -p 8180:80 -p 8121:21 -p 8122:22 --rm -v ${NGD_BALDER_BASE}/volumes/galaxy_storage/:/export  bgruening/galaxy-stable
  elif [ "$(whoami)" == "rnebot" ] ; then
    ssh rnebot@balder docker run --name galaxy_devel_rnebot -d -p 8180:80 -p 8121:21 -p 8122:22 --rm -v /home/rnebot/DATOS/galaxy_storage/:/export  bgruening/galaxy-stable
  elif [ "$(whoami)" == "acurbelo" ] ; then
#    ssh acurbelo@balder docker run --name galaxy_devel_acurbelo -d -p 8280:80 -p 8221:21 -p 8222:22 --rm -v /var/lib/nextgendem/galaxy_storage/:/export  bgruening/galaxy-stable
    docker run --name galaxy_devel -d -p 8080:80 -p 8021:21 -p 8022:22 --rm -v /var/lib/nextgendem/galaxy_storage/:/export  bgruening/galaxy-stable:20.05
  elif [ "$(whoami)" == "paula" ] ; then
    docker run --name galaxy_devel -d -p 8080:80 -p 8021:21 -p 8022:22 --rm -v /home/paula/galaxy_storage/:/export  bgruening/galaxy-stable
  elif [ "$(whoami)" == "daniel" ] ; then
    #ssh dreyes@balder docker run --name galaxy_devel_dreyes -d -p 8480:80 -p 8421:21 -p 8422:22 --rm -v /home/daniel/Documentos/DATOS/galaxy_storage/:/export bgruening/galaxy-stable
    #docker run --name galaxy_devel -d -p 8080:80 --privileged=true -v /home/daniel/Documentos/galaxy_storage/:/export/ bgruening/galaxy-stable:latest
    docker run --name galaxy_devel -d -p 8080:80 -p 8021:21 -p 8022:22 --rm -v /home/daniel/Documentos/galaxy_storage/:/export  bgruening/galaxy-stable:latest
  fi
fi

# postgis
postgis_started="yes"
if [ "$(hostname)" == "balder" ] ; then
  postgis_started=$(docker ps -q -f name="^postgis$")
elif [ "$(whoami)" == "rnebot" ] ; then
  postgis_started=$(docker ps -q -f name=postgis)
elif [ "$(whoami)" == "acurbelo" ] ; then
  postgis_started=$(docker ps -q -f name=postgis)
elif [ "$(whoami)" == "paula" ] ; then
  postgis_started=$(docker ps -q -f name=postgis)
elif [ "$(whoami)" == "daniel" ] ; then
  postgis_started=$(docker ps -q -f name=postgis)
fi

if [ ! "$postgis_started" ] ; then
  echo Starting Postgis
  if [ "$(hostname)" == "balder" ] ; then
    docker run -d -v ${NGD_BALDER_BASE}/volumes/geoserver_pg:/var/lib/postgresql -p 7434:5432 --rm -e POSTGRES_DBNAME=ngd_geoserver -e POSTGRES_PASS=postgres -e POSTGRES_USER=postgres -e ALLOW_IP_RANGE=0.0.0.0/0 --name postgis -t kartoza/postgis:13.0
  elif [ "$(whoami)" == "rnebot" ] ; then
    docker run -d -v /home/rnebot/DATOS/geoserver/pg:/var/lib/postgresql -p 5435:5432 --rm -e POSTGRES_DBNAME=ngd_geoserver -e POSTGRES_PASS=postgres -e POSTGRES_USER=postgres -e ALLOW_IP_RANGE=0.0.0.0/0 --name postgis -t kartoza/postgis:13.0
  elif [ "$(whoami)" == "acurbelo" ] ; then
    docker run -d -v /var/lib/nextgendem/geoserver/pg:/var/lib/postgresql -p 5435:5432 --rm -e POSTGRES_DBNAME=ngd_geoserver -e POSTGRES_PASS=postgres -e POSTGRES_USER=postgres -e ALLOW_IP_RANGE=0.0.0.0/0 --name postgis -t kartoza/postgis:13.0
  elif [ "$(whoami)" == "paula" ] ; then
    docker run -d -v /home/paula/DATOS/geoserver/pg:/var/lib/postgresql -p 5435:5432 --rm -e POSTGRES_DBNAME=ngd_geoserver -e POSTGRES_PASS=postgres -e POSTGRES_USER=postgres -e ALLOW_IP_RANGE=0.0.0.0/0 --name postgis -t kartoza/postgis:13.0
  elif [ "$(whoami)" == "daniel" ] ; then
    docker run -d -v /home/daniel/Documentos/DATOS/geoserver/pg:/var/lib/postgresql -p 5435:5432 --rm -e POSTGRES_DBNAME=ngd_geoserver -e POSTGRES_PASS=postgres -e POSTGRES_USER=postgres -e ALLOW_IP_RANGE=0.0.0.0/0 --name postgis -t kartoza/postgis:13.0
  elif [ "$(whoami)" == "carlos" ] ; then
    docker run -d -v /home/carlos/DATOS/bcs/pg:/var/lib/postgresql -p 5435:5432 --rm -e POSTGRES_DBNAME=ngd_geoserver -e POSTGRES_PASS=postgres -e POSTGRES_USER=postgres -e ALLOW_IP_RANGE=0.0.0.0/0 --name postgis -t kartoza/postgis:13.0
  fi
fi


if [ "$(hostname)" == "balder" ] ; then
  if [ ! -d "${NGD_BALDER_BASE}/volumes/geoserver" ]; then
    mkdir ${NGD_BALDER_BASE}/volumes/geoserver/; sudo chown 1000 ${NGD_BALDER_BASE}/volumes/geoserver
  fi
elif [ "$(whoami)" == "rnebot" ] ; then
  if [ ! -d "/home/rnebot/DATOS/geoserver/geoserver" ]; then
    mkdir /home/rnebot/DATOS/geoserver/geoserver/; sudo chown 1000 /home/rnebot/DATOS/geoserver/geoserver
  fi
elif [ "$(whoami)" == "acurbelo" ] ; then
  if [ ! -d "/var/lib/nextgendem/geoserver/geoserver" ]; then
    mkdir /var/lib/nextgendem/geoserver/geoserver; sudo chown 1000 /var/lib/nextgendem/geoserver/geoserver
  fi
elif [ "$(whoami)" == "paula" ] ; then
   if [ ! -d "/home/paula/DATOS/geoserver/geoserver" ]; then
    mkdir /home/paula/DATOS/geoserver/geoserver; sudo chown 1000 /home/paula/DATOS/geoserver/geoserver
   fi
elif [ "$(whoami)" == "daniel" ] ; then
  if [ ! -d "/home/daniel/Documentos/DATOS/geoserver/geoserver" ]; then
    mkdir /home/daniel/Documentos/DATOS/geoserver/geoserver; sudo chown 1000 /home/daniel/Documentos/DATOS/geoserver/geoserver
  fi
elif [ "$(whoami)" == "carlos" ] ; then
  if [ ! -d "/home/carlos/DATOS/bcs/geoserver" ]; then
    mkdir /home/carlos/DATOS/bcs/geoserver; sudo chown 1000 /home/carlos/DATOS/bcs/geoserver
  fi
fi


# Geoserver
geoserver_started="yes"
if [ "$(hostname)" == "balder" ] ; then
    geoserver_started=$(docker ps -q -f name="^geoserver$")
elif [ "$(whoami)" == "rnebot" ] ; then
    geoserver_started=$(docker ps -q -f name=geoserver)
elif [ "$(whoami)" == "acurbelo" ] ; then
    geoserver_started=$(docker ps -q -f name=geoserver)
elif [ "$(whoami)" == "paula" ] ; then
    geoserver_started=$(docker ps -q -f name=geoserver)
elif [ "$(whoami)" == "daniel" ] ; then
    geoserver_started=$(docker ps -q -f name=geoserver)
elif [ "$(whoami)" == "carlos" ] ; then
    geoserver_started=$(docker ps -q -f name=geoserver)
fi

if [ ! "$geoserver_started" ] ; then
  echo Starting Geoserver
  if [ "$(hostname)" == "balder" ] ; then
    docker run --name geoserver --link postgis:postgis -p 8280:8080 --rm -v ${NGD_BALDER_BASE}/volumes/geoserver/:/opt/geoserver/data_dir -d -e DB_BACKEND=POSTGRES -e HOST=postgis -e POSTGRES_PORT=5432 -e POSTGRES_DB=geoserver -e POSTGRES_USER=postgres -e POSTGRES_PASS=postgres -e USERNAME=postgres -e PASS=postgres -e GEOSERVER_ADMIN_PASSWORD=ngd_ad37 -e GEOSERVER_ADMIN_USER=admin kartoza/geoserver:2.19.0
  elif [ "$(whoami)" == "rnebot" ] ; then
    docker run --name "geoserver"   --link postgis:postgis  -p 9180:8080 --rm -v /home/rnebot/DATOS/geoserver/geoserver/:/opt/geoserver/data_dir  -d  -e DB_BACKEND=POSTGRES  -e HOST=postgis  -e POSTGRES_PORT=5432  -e POSTGRES_DB=geoserver  -e POSTGRES_USER=postgres  -e POSTGRES_PASS=postgres  -e USERNAME=postgres  -e PASS=postgres  -e GEOSERVER_ADMIN_PASSWORD=ngd_ad37   -e GEOSERVER_ADMIN_USER=admin  kartoza/geoserver:2.19.0
  elif [ "$(whoami)" == "acurbelo" ] ; then
    docker run --name "geoserver"  --link postgis:postgis  -p 9180:8080 --rm -v /var/lib/nextgendem/geoserver/geoserver:/opt/geoserver/data_dir  -d  -e DB_BACKEND=POSTGRES  -e HOST=postgis  -e POSTGRES_PORT=5432  -e POSTGRES_DB=geoserver  -e POSTGRES_USER=postgres  -e POSTGRES_PASS=postgres  -e USERNAME=postgres  -e PASS=postgres  -e GEOSERVER_ADMIN_PASSWORD=ngd_ad37   -e GEOSERVER_ADMIN_USER=admin  kartoza/geoserver:2.19.0
  elif [ "$(whoami)" == "paula" ] ; then
    docker run --name "geoserver"  --link postgis:postgis  -p 9180:8080 --rm -v /home/paula/DATOS/geoserver/geoserver:/opt/geoserver/data_dir  -d  -e DB_BACKEND=POSTGRES  -e HOST=postgis  -e POSTGRES_PORT=5432  -e POSTGRES_DB=ngd_eoserver  -e POSTGRES_USER=postgres  -e POSTGRES_PASS=postgres  -e USERNAME=postgres  -e PASS=postgres  -e GEOSERVER_ADMIN_PASSWORD=ngd_ad37   -e GEOSERVER_ADMIN_USER=admin  kartoza/geoserver:2.19.0
  elif [ "$(whoami)" == "daniel" ] ; then
    docker run --name "geoserver"  --link postgis:postgis  -p 9180:8080 --rm -v /home/daniel/Documentos/DATOS/geoserver/geoserver:/opt/geoserver/data_dir  -d  -e DB_BACKEND=POSTGRES  -e HOST=postgis  -e POSTGRES_PORT=5432  -e POSTGRES_DB=geoserver  -e POSTGRES_USER=postgres  -e POSTGRES_PASS=postgres  -e USERNAME=postgres  -e PASS=postgres  -e GEOSERVER_ADMIN_PASSWORD=ngd_ad37   -e GEOSERVER_ADMIN_USER=admin  kartoza/geoserver:2.19.0
  elif [ "$(whoami)" == "carlos" ] ; then
    docker run --name "geoserver"  --link postgis:postgis  -p 9180:8080 --rm -v /home/carlos/DATOS/bcs/geoserver:/opt/geoserver/data_dir  -d  -e DB_BACKEND=POSTGRES  -e HOST=postgis  -e POSTGRES_PORT=5432  -e POSTGRES_DB=geoserver  -e POSTGRES_USER=postgres  -e POSTGRES_PASS=postgres  -e USERNAME=postgres  -e PASS=postgres  -e GEOSERVER_ADMIN_PASSWORD=ngd_ad37   -e GEOSERVER_ADMIN_USER=admin  kartoza/geoserver:2.19.0
  fi
fi

# CD to BCS-BACKEND source code (needed for proper Celery execution)
CWD=$(basename "$PWD")
if test $CWD = "docker_init"
then
  cd ..
fi

# BACKEND IMAGE (Backend+Celery) or just Celery (for local development)
if [ "$(hostname)" == "balder" ] ; then
  # --rm
  balder_started=$(docker ps -q -f name="^ngd_backend$")
  if [ ! "$balder_started" ] ; then
    echo Starting Backend
    docker run --name ngd_backend -p 80:80 --rm \
        --add-host balder.itccanarias.org:172.17.0.1 \
        -e BCS_CONFIG_FILE="/root/.local/share/ngd-backend/ngd_local.conf" \
        --link postgres_devel_2:app_db --link postgres_devel:chado_db --link postgis:postgis \
        --link redis:redis --link geoserver:geoserver \
        -v ${NGD_BALDER_BASE}/private-conf/bcs_balder_dev.conf:/root/.local/share/ngd-backend/ngd_local.conf \
        -v ${NGD_BALDER_BASE}/private-conf/firebase-key.json:/root/.local/share/ngd-backend/firebase-key.json \
        -v ${NGD_BALDER_BASE}/private-conf/resources.json:/root/.local/share/ngd-backend/resources.json \
        -v ${NGD_BALDER_BASE}/private-conf/.ssh/:/root/.ssh/ \
        -d nextgendem-mac/ngd-bcs-backend:latest
  fi
else
  # Celery
  if [ "$(whoami)" == "rnebot" ] ; then
    if [ "$#" -gt 0 ] ; then
      start_celery="yes"
    else
      start_celery="no"
    fi
  else
    start_celery="yes"
  fi

  if [ "$start_celery" == "yes" ] ; then
    # Worker AND Beat (only for development; NOT for production -use Supervisor and two separate processes-)
    celery -A biobarcoding.tasks worker --beat --loglevel=info
  fi
fi
