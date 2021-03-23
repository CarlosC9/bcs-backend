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
    #init_chado "$(dirname $0)"
  fi
fi

# Galaxy
# api key = fakekey; user = admin; password = password
galaxy_started="yes"
if [ "$(whoami)" == "rnebot" ] && [ "$#" -gt 0 ] ; then
  galaxy_started=$(ssh rnebot@balder docker ps -q -f name=galaxy_devel_rnebot)
elif [ "$(whoami)" == "acurbelo" ] ; then
  galaxy_started=$(ssh acurbelo@balder docker ps -q -f name=galaxy_devel_acurbelo)
#  galaxy_started=$(docker ps -q -f name=galaxy_devel)
elif [ "$(whoami)" == "paula" ] ; then
  galaxy_started=$(docker ps -q -f name=galaxy_devel)
elif [ "$(whoami)" == "daniel" ] ; then
  galaxy_started=$(docker ps -q -f name=galaxy_devel)
  #galaxy_started=$(ssh dreyes@balder docker ps -q -f name=galaxy_devel_dreyes)
fi

if [ ! $galaxy_started ] ; then
  echo Starting Galaxy
  if [ "$(whoami)" == "rnebot" ] ; then
    ssh rnebot@balder docker run --name galaxy_devel_rnebot -d -p 8180:80 -p 8121:21 -p 8122:22 --rm -v /home/rnebot/DATOS/galaxy_storage/:/export  bgruening/galaxy-stable
  elif [ "$(whoami)" == "acurbelo" ] ; then
    ssh acurbelo@balder docker run --name galaxy_devel_acurbelo -d -p 8280:80 -p 8221:21 -p 8222:22 --rm -v /var/lib/nextgendem/galaxy_storage/:/export  bgruening/galaxy-stable
#    docker run --name galaxy_devel -d -p 8080:80 -p 8021:21 -p 8022:22 --rm -v /var/lib/nextgendem/galaxy_storage/:/export  bgruening/galaxy-stable
  elif [ "$(whoami)" == "paula" ] ; then
    docker run --name galaxy_devel -d -p 8080:80 -p 8021:21 -p 8022:22 --rm -v /home/paula/galaxy_storage/:/export  bgruening/galaxy-stable
  elif [ "$(whoami)" == "daniel" ] ; then
    #ssh dreyes@balder docker run --name galaxy_devel_dreyes -d -p 8480:80 -p 8421:21 -p 8422:22 --rm -v /home/daniel/Documentos/DATOS/galaxy_storage/:/export bgruening/galaxy-stable
    #docker run --name galaxy_devel -d -p 8080:80 --privileged=true -v /home/daniel/Documentos/galaxy_storage/:/export/ bgruening/galaxy-stable:latest
    docker run --name galaxy_devel -d -p 8080:80 -p 8021:21 -p 8022:22 --rm -v /home/daniel/Documentos/galaxy_storage/:/export  bgruening/galaxy-stable:latest
  fi
fi

# Geoserver
geoserver_started="yes"
if [ "$(whoami)" == "rnebot" ] && [ "$#" -gt 0 ] ; then
  geoserver_started=$(ssh rnebot@balder docker ps -q -f name=geoserver_devel_rnebot)
elif [ "$(whoami)" == "acurbelo" ] ; then
  geoserver_started=$(docker ps -q -f name=geoserver_devel)
elif [ "$(whoami)" == "paula" ] ; then
  geoserver_started=$(docker ps -q -f name=geoserver_devel)
elif [ "$(whoami)" == "daniel" ] ; then
  geoserver_started=$(ssh dreyes@balder docker ps -q -f name=geoserver_devel_dreyes)
fi

if [ ! geoserver_started ] ; then
  echo Starting Geoserver
  if [ "$(whoami)" == "rnebot" ] ; then
    ssh rnebot@balder docker run --name geoserver_devel_rnebot -d -p 9180:80 --rm -v /home/rnebot/DATOS/geoserver_storage/:/export  bgruening/galaxy-stable
  elif [ "$(whoami)" == "acurbelo" ] ; then
    docker run --name galaxy_devel -d -p 8080:80 -p 8021:21 -p 8022:22 --rm -v /var/lib/nextgendem/galaxy_storage/:/export  bgruening/galaxy-stable
  elif [ "$(whoami)" == "paula" ] ; then
    docker run --name galaxy_devel -d -p 8080:80 -p 8021:21 -p 8022:22 --rm -v /home/paula/galaxy_storage/:/export  bgruening/galaxy-stable
  elif [ "$(whoami)" == "daniel" ] ; then
    #docker run --name galaxy_devel -d -p 8080:80 -p 8021:21 -p 8022:22 --rm -v /home/daniel/Documentos/DATOS/galaxy_storage/:/export  bgruening/galaxy-stable
    ssh dreyes@balder docker run --name galaxy_devel_dreyes -d -p 8480:80 -p 8421:21 -p 8422:22 --rm -v /home/daniel/Documentos/galaxy_data/:/export bgruening/galaxy-stable
  #docker run -d -p 8080:80 -p 8021:21 -p 8800:8800 --privileged=true -v /home/daniel/Documentos/galaxy_storage/:/export/ bgruening/galaxy-stable
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
