#!/bin/bash
# REDIS
if [ ! "$(docker ps -q -f name=redis)" ] ; then
  echo Starting REDIS
  docker run --name redis --rm -d -p 6379:6379 redis
fi

# Initialize EDAM ontology, NCBI taxonomy, ...
function init_chado {
  echo Waiting for PostgreSQL-Chado
  pg_isready -d postgres -h localhost -p 5432 -U postgres
  while [ $? -ne 0 ]
  do
    echo PostgreSQL-Chado is not ready yet.
    sleep 3
    pg_isready -d postgres -h localhost -p 5432 -U postgres
  done
  echo Initializing Chado
  cd $1/docker_init
  ./init.sh
}


if [ "$(whoami)" == "rnebot" ] && [ "$#" -gt 0 ] ; then
  echo "TODO: INICIALIZAR VARIABLES DE ENTORNO!"
elif [ "$(whoami)" == "acurbelo" ] ; then
  echo "TODO: INICIALIZAR VARIABLES DE ENTORNO!"
elif [ "$(whoami)" == "paula" ] ; then
  export ENDPOINT_URL="http://localhost:5000"
  export COOKIES_FILE_PATH="/home/paula/Documentos/NEXTGENDEM/curl/bcs-cookies.txt"
elif [ "$(whoami)" == "daniel" ] ; then
  export ENDPOINT_URL="http://localhost:5000"
  export COOKIES_FILE_PATH="/home/daniel/Documentos/Projects/curl/bcs-cookies.txt"
fi

# PostgreSQL
if [ ! "$(docker ps -q -f name=postgres_devel)" ] ; then
  echo Starting PostgreSQL-Chado
  if [ "$(whoami)" == "rnebot" ] ; then
    docker run --name postgres_devel -d -p 5432:5432 --rm -e POSTGRES_PASSWORD=postgres -e INSTALL_CHADO_SCHEMA=1 -e INSTALL_YEAST_DATA=0 -e PGDATA=/var/lib/postgresql/data/ -v /home/rnebot/DATOS/pg_devel:/var/lib/postgresql/data quay.io/galaxy-genome-annotation/chado:1.31-jenkins97-pg9.5
    init_chado "$(dirname $0)"
#    docker run --name postgres_devel -d -p 5432:5432 --rm -e POSTGRES_PASSWORD=postgres -v /home/rnebot/DATOS/pg_devel:/var/lib/postgresql/data postgres
  elif [ "$(whoami)" == "acurbelo" ] ; then
    docker run --name postgres_devel -d -p 5432:5432 --rm -e POSTGRES_PASSWORD=postgres -e INSTALL_CHADO_SCHEMA=1 -e INSTALL_YEAST_DATA=0 -e PGDATA=/var/lib/postgresql/data/ -v /var/lib/nextgendem/pg_devel:/var/lib/postgresql/data quay.io/galaxy-genome-annotation/chado:1.31-jenkins97-pg9.5
    init_chado "$(dirname $0)"
  elif [ "$(whoami)" == "paula" ] ; then
    docker run --name postgres_devel -d -p 5432:5432 --rm -e POSTGRES_PASSWORD=postgres -e INSTALL_CHADO_SCHEMA=1 -e INSTALL_YEAST_DATA=0 -e PGDATA/var/lib/postgresql/data/ -v /home/paula/DATOS/pg_devel:/var/lib/postgresql/data quay.io/galaxy-genome-annotation/chado:1.31-jenkins97-pg9.5
  elif [ "$(whoami)" == "daniel" ] ; then
    docker run --name postgres_devel -d -p 5432:5432 --rm -e POSTGRES_PASSWORD=postgres -e INSTALL_CHADO_SCHEMA=1 -e INSTALL_YEAST_DATA=0 -e PGDATA=/var/lib/postgresql/data/ -v /home/daniel/Documentos/DATOS/pg_devel:/var/lib/postgresql/data quay.io/galaxy-genome-annotation/chado:1.31-jenkins97-pg9.5
    init_chado "$(dirname $0)"
  fi
fi

# Galaxy
# api key = fakekey; user = admin; password = password
galaxy_started="yes"
if [ "$(whoami)" == "rnebot" ] && [ "$#" -gt 0 ] ; then
  galaxy_started=$(ssh rnebot@balder docker ps -q -f name=galaxy_devel_rnebot)
elif [ "$(whoami)" == "acurbelo" ] ; then
  galaxy_started=$(ssh acurbelo@balder docker ps -q -f name=galaxy_devel_acurbelo)
elif [ "$(whoami)" == "paula" ] ; then
  galaxy_started=$(docker ps -q -f name=galaxy_devel)
elif [ "$(whoami)" == "daniel" ] ; then
  galaxy_started=$(docker ps -q -f name=galaxy_devel)
  #galaxy_started=$(ssh dreyes@balder docker ps -q -f name=galaxy_devel_dreyes)
fi

if [ ! "$galaxy_started" ] ; then
  echo Starting Galaxy
  if [ "$(whoami)" == "rnebot" ] ; then
    ssh rnebot@balder docker run --name galaxy_devel_rnebot -d -p 8180:80 -p 8121:21 -p 8122:22 --rm -v /home/rnebot/DATOS/galaxy_storage/:/export  bgruening/galaxy-stable
  elif [ "$(whoami)" == "acurbelo" ] ; then
    ssh acurbelo@balder docker run --name galaxy_devel_acurbelo -d -p 8280:80 -p 8221:21 -p 8222:22 --rm -v /var/lib/nextgendem/galaxy_storage/:/export  bgruening/galaxy-stable
  elif [ "$(whoami)" == "paula" ] ; then
    docker run --name galaxy_devel -d -p 8080:80 -p 8021:21 -p 8022:22 --rm -v /home/paula/galaxy_storage/:/export  bgruening/galaxy-stable
  elif [ "$(whoami)" == "daniel" ] ; then
    #ssh dreyes@balder docker run --name galaxy_devel_dreyes -d -p 8480:80 -p 8421:21 -p 8422:22 --rm -v /home/daniel/Documentos/DATOS/galaxy_storage/:/export bgruening/galaxy-stable
    #docker run --name galaxy_devel -d -p 8080:80 --privileged=true -v /home/daniel/Documentos/galaxy_storage/:/export/ bgruening/galaxy-stable:latest
    docker run --name galaxy_devel -d -p 8080:80 -p 8021:21 -p 8022:22 --rm -v /home/daniel/Documentos/galaxy_storage/:/export  bgruening/galaxy-stable:latest
  fi
fi

if [ "$(whoami)" == "rnebot" ] ; then
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
  if [ ! -d "daniel/Documentos/geoserver" ]; then
    mkdir /home/daniel/Documentos/geoserver
  fi
fi

# postgis
postgis_started="yes"
if [ "$(whoami)" == "rnebot" ] && [ "$#" -gt 0 ] ; then
  postgis_started=$(docker ps -q -f name=postgis)
elif [ "$(whoami)" == "acurbelo" ] ; then
  postgis_started=$(docker ps -q -f name=postgis)
elif [ "$(whoami)" == "paula" ] ; then
  postgis_started=$(docker ps -q -f name=postgis)
elif [ "$(whoami)" == "daniel" ] ; then
  postgis_started=$(docker ps -q -f name=postgis)
fi

if [ ! $postgis_started ] ; then
  echo Starting Postgis
  if [ "$(whoami)" == "rnebot" ] ; then
    docker run -d -v /home/rnebot/DATOS/geoserver/pg:/var/lib/postgresql -p 5435:5432 --rm -e POSTGRES_DBNAME=ngd_geoserver -e POSTGRES_PASS=postgres -e POSTGRES_USER=postgres -e ALLOW_IP_RANGE=0.0.0.0/0 --name postgis -t kartoza/postgis:13.0
  elif [ "$(whoami)" == "acurbelo" ] ; then
    docker run -d -v /var/lib/nextgendem/geoserver/pg:/var/lib/postgresql -p 5435:5432 --rm -e POSTGRES_DBNAME=ngd_geoserver -e POSTGRES_PASS=postgres -e POSTGRES_USER=postgres -e ALLOW_IP_RANGE=0.0.0.0/0 --name postgis -t kartoza/postgis:13.0
  elif [ "$(whoami)" == "paula" ] ; then
    docker run -d -v /home/paula/DATOS/geoserver/pg:/var/lib/postgresql -p 5435:5432 --rm -e POSTGRES_DBNAME=ngd_geoserver -e POSTGRES_PASS=postgres -e POSTGRES_USER=postgres -e ALLOW_IP_RANGE=0.0.0.0/0 --name postgis -t kartoza/postgis:13.0
  elif [ "$(whoami)" == "daniel" ] ; then
    docker run -d -v /home/daniel/Documentos/geoserver/pg:/var/lib/postgresql -p 5435:5432 --rm -e POSTGRES_DBNAME=ngd_geoserver -e POSTGRES_PASS=postgres -e POSTGRES_USER=postgres -e ALLOW_IP_RANGE=0.0.0.0/0 --name postgis -t kartoza/postgis:13.0
  fi
fi


if [ "$(whoami)" == "rnebot" ] ; then
  if [ ! -d "/home/rnebot/DATOS/geoserver/geoserver" ]; then
    mkdir /home/rnebot/DATOS/geoserver/geoserver/ && sudo chown 1000 /home/rnebot/DATOS/geoserver/geoserver
  fi
elif [ "$(whoami)" == "acurbelo" ] ; then
  if [ ! -d "/var/lib/nextgendem/geoserver/geoserver" ]; then
    mkdir /var/lib/nextgendem/geoserver/geoserver && sudo chown 1000 /var/lib/nextgendem/geoserver/geoserver
  fi
elif [ "$(whoami)" == "paula" ] ; then
   if [ ! -d "/home/paula/DATOS/geoserver/geoserver" ]; then
    mkdir /home/paula/DATOS/geoserver/geoserver && sudo chown 1000 /home/paula/DATOS/geoserver/geoserver
   fi
elif [ "$(whoami)" == "daniel" ] ; then
  if [ ! -d "daniel/Documentos/geoserver/geoserver" ]; then
    mkdir /home/daniel/Documentos/geoserver/geoserver && sudo chown 1000 /home/rnebot/DATOS/geoserver/geoserver
  fi
fi


# Geoserver
geoserver_started="yes"
if [ "$(whoami)" == "rnebot" ] && [ "$#" -gt 0 ] ; then
    geoserver_started=$(docker ps -q -f name=geoserver)
  elif [ "$(whoami)" == "acurbelo" ] ; then
    geoserver_started=$(docker ps -q -f name=geoserver)
  elif [ "$(whoami)" == "paula" ] ; then
    geoserver_started=$(docker ps -q -f name=geoserver)
  elif [ "$(whoami)" == "daniel" ] ; then
    geoserver_started=$(docker ps -q -f name=geoserver)
fi

if [ ! $geoserver_started ] ; then
  echo Starting Geoserver
  if [ "$(whoami)" == "rnebot" ] ; then
    docker run --name "geoserver"   --link postgis:postgis  -p 9180:8080 --rm -v /home/rnebot/DATOS/geoserver/geoserver/:/opt/geoserver/data_dir  -d  -e DB_BACKEND=POSTGRES  -e HOST=postgis  -e POSTGRES_PORT=5432  -e POSTGRES_DB=geoserver  -e POSTGRES_USER=postgres  -e POSTGRES_PASS=postgres  -e USERNAME=postgres  -e PASS=postgres  -e GEOSERVER_ADMIN_PASSWORD=ngd_ad37   -e GEOSERVER_ADMIN_USER=ngd_admin  kartoza/geoserver:2.19.0
  elif [ "$(whoami)" == "acurbelo" ] ; then
    docker run --name "geoserver"  --link postgis:postgis  -p 9180:8080 --rm -v /var/lib/nextgendem/geoserver/geoserver:/opt/geoserver/data_dir  -d  -e DB_BACKEND=POSTGRES  -e HOST=postgis  -e POSTGRES_PORT=5432  -e POSTGRES_DB=geoserver  -e POSTGRES_USER=postgres  -e POSTGRES_PASS=postgres  -e USERNAME=postgres  -e PASS=postgres  -e GEOSERVER_ADMIN_PASSWORD=ngd_ad37   -e GEOSERVER_ADMIN_USER=ngd_admin  kartoza/geoserver:2.19.0
  elif [ "$(whoami)" == "paula" ] ; then
    docker run --name "geoserver"  --link postgis:postgis  -p 9180:8080 --rm -v /home/paula/DATOS/geoserver/geoserver:/opt/geoserver/data_dir  -d  -e DB_BACKEND=POSTGRES  -e HOST=postgis  -e POSTGRES_PORT=5432  -e POSTGRES_DB=ngd_eoserver  -e POSTGRES_USER=postgres  -e POSTGRES_PASS=postgres  -e USERNAME=postgres  -e PASS=postgres  -e GEOSERVER_ADMIN_PASSWORD=ngd_ad37   -e GEOSERVER_ADMIN_USER=ngd_admin  kartoza/geoserver:2.19.0
  elif [ "$(whoami)" == "daniel" ] ; then
    docker run --name "geoserver"  --link postgis:postgis  -p 9180:8080 --rm -v /home/daniel/Documentos/geoserver/geoserver:/opt/geoserver/data_dir  -d  -e DB_BACKEND=POSTGRES  -e HOST=postgis  -e POSTGRES_PORT=5432  -e POSTGRES_DB=geoserver  -e POSTGRES_USER=postgres  -e POSTGRES_PASS=postgres  -e USERNAME=postgres  -e PASS=postgres  -e GEOSERVER_ADMIN_PASSWORD=ngd_ad37   -e GEOSERVER_ADMIN_USER=ngd_admin  kartoza/geoserver:2.19.0
  fi
fi

# CD to BCS-BACKEND source code (needed for proper Celery execution)
CWD=$(basename "$PWD")
if test $CWD = "docker_init"
then
  cd ..
fi

# Worker AND Beat (only for development; NOT for production -use Supervisor and two separate processes-)
celery -A biobarcoding.tasks worker --beat --loglevel=info
