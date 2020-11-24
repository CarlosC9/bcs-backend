#!/bin/bash
# REDIS
if [ ! "$(docker ps -q -f name=redis)" ] ; then
  echo Starting REDIS
  docker run --name redis --rm -d -p 6379:6379 redis
fi

# PostgreSQL
if [ ! "$(docker ps -q -f name=postgres_devel)" ] ; then
  echo Starting PostgreSQL-Chado
  if [ "$(whoami)" == "rnebot" ] ; then
    docker run --name postgres_devel -d -p 5432:5432 --rm -e POSTGRES_PASSWORD=postgres -e INSTALL_CHADO_SCHEMA=1 -e INSTALL_YEAST_DATA=0 -e PGDATA=/var/lib/postgresql/data/ -v /home/rnebot/DATOS/pg_devel:/var/lib/postgresql/data quay.io/galaxy-genome-annotation/chado:1.31-jenkins97-pg9.5
#    docker run --name postgres_devel -d -p 5432:5432 --rm -e POSTGRES_PASSWORD=postgres -v /home/rnebot/DATOS/pg_devel:/var/lib/postgresql/data postgres
  elif [ "$(whoami)" == "acurbelo" ] ; then
    docker run --name postgres_devel -d -p 5432:5432 --rm -e POSTGRES_PASSWORD=postgres -e INSTALL_CHADO_SCHEMA=1 -e INSTALL_YEAST_DATA=0 -e PGDATA=/var/lib/postgresql/data/ -v /var/lib/nextgendem/pg_devel:/var/lib/postgresql/data quay.io/galaxy-genome-annotation/chado:1.31-jenkins97-pg9.5
  elif [ "$(whoami)" == "paula" ] ; then
    docker run --name postgres_devel -d -p 5432:5432 --rm -e POSTGRES_PASSWORD=postgres -e INSTALL_CHADO_SCHEMA=1 -e INSTALL_YEAST_DATA=0 -e PGDATA/var/lib/postgresql/data/ -v /home/paula/DATOS/pg_devel:/var/lib/postgresql/data quay.io/galaxy-genome-annotation/chado:1.31-jenkins97-pg9.5
  elif [ "$(whoami)" == "daniel" ] ; then
    docker run --name postgres_devel -d -p 5432:5432 --rm -e POSTGRES_PASSWORD=postgres -e INSTALL_CHADO_SCHEMA=1 -e INSTALL_YEAST_DATA=0 -e PGDATA=/var/lib/postgresql/data/ -v /home/daniel/Documentos/DATOS/pg_devel:/var/lib/postgresql/data quay.io/galaxy-genome-annotation/chado:1.31-jenkins97-pg9.5
    timeout 30 ping 8.8.8.8 > NUL #we need to wait 30s until the database is created
    cd ~/Documentos/GIT/bcs-backend/docker_init
    ./init.sh
  fi
fi



# Galaxy
# api key = fakekey; user = admin; password = password
galaxy_started="yes"
if [ "$(whoami)" == "rnebot" ] && [ "$#" -gt 0 ] ; then
  galaxy_started=$(ssh rnebot@balder docker ps -q -f name=galaxy_devel_rnebot)
elif [ "$(whoami)" == "acurbelo" ] ; then
  galaxy_started=$(docker ps -q -f name=galaxy_devel)
elif [ "$(whoami)" == "paula" ] ; then
  galaxy_started=$(docker ps -q -f name=galaxy_devel)
elif [ "$(whoami)" == "daniel" ] ; then
  galaxy_started=$(ssh dreyes@balder docker ps -q -f name=galaxy_devel_dreyes)
fi

if [ ! $galaxy_started ] ; then
  echo Starting Galaxy
  if [ "$(whoami)" == "rnebot" ] ; then
    ssh rnebot@balder docker run --name galaxy_devel_rnebot -d -p 8180:80 -p 8121:21 -p 8122:22 --rm -v /home/rnebot/DATOS/galaxy_storage/:/export  bgruening/galaxy-stable
  elif [ "$(whoami)" == "acurbelo" ] ; then
    docker run --name galaxy_devel -d -p 8080:80 -p 8021:21 -p 8022:22 --rm -v /var/lib/nextgendem/galaxy_storage/:/export  bgruening/galaxy-stable
  elif [ "$(whoami)" == "paula" ] ; then
    docker run --name galaxy_devel -d -p 8080:80 -p 8021:21 -p 8022:22 --rm -v /home/paula/galaxy_storage/:/export  bgruening/galaxy-stable
  elif [ "$(whoami)" == "daniel" ] ; then
    ssh dreyes@balder docker run --name galaxy_devel_dreyes -d -p 8480:80 -p 8421:21 -p 8422:22 --rm -v /home/daniel/Documentos/DATOS/galaxy_storage/:/export bgruening/galaxy-stable
    #echo El galaxy no funciona
  fi
fi

# CD to BCS-BACKEND source code (needed for proper Celery execution)
if [ "$(whoami)" == "rnebot" ] ; then
  cd ~/GoogleDrive/AA_NEXTGENDEM/bcs-backend/
elif [ "$(whoami)" == "acurbelo" ] ; then
  cd ~/Proyectos/NEXTGENDEM/bcs-backend/
elif [ "$(whoami)" == "paula" ] ; then
  cd ~/Documentos/NEXTGENDEM/bcs/bcs-backend/
elif [ "$(whoami)" == "daniel" ] ; then
  cd ~/Documentos/GIT/bcs-backend/
fi

# CELERY
# Worker AND Beat (only for development; NOT for production -use Supervisor and two separate processes-)
celery -A biobarcoding.tasks.definitions worker --beat --loglevel=info
