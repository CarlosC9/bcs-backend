#!/bin/bash
# REDIS
if [ ! "$(docker ps -q -f name=redis)" ] ; then
  docker run --name redis --rm -d -p 6379:6379 redis
fi

# PostgreSQL
if [ ! "$(docker ps -q -f name=postgres_devel)" ] ; then
  if [ "$(whoami)" == "rnebot" ] ; then
    docker run --name postgres_devel -d -p 5432:5432 --rm -e POSTGRES_PASSWORD=postgres -v /home/rnebot/DATOS/pg_devel:/var/lib/postgresql/data postgres
  elif [ "$(whoami)" == "acurbelo"] ; then
    docker run --name postgres_devel -d -p 5432:5432 --rm -e POSTGRES_PASSWORD=postgres -v <........>:/var/lib/postgresql/data postgres
  elif [ "$(whoami)" == "pmoreno" ] ; then
    docker run --name postgres_devel -d -p 5432:5432 --rm -e POSTGRES_PASSWORD=postgres -v <........>:/var/lib/postgresql/data postgres
  elif [ "$(whoami)" == "dreyes" ] ; then
    docker run --name postgres_devel -d -p 5432:5432 --rm -e POSTGRES_PASSWORD=postgres -v <........>:/var/lib/postgresql/data postgres
  fi
fi

# CD to BCS-BACKEND source code (needed for proper Celery execution)
if [ "$(whoami)" == "rnebot" ] ; then
  cd ~/GoogleDrive/AA_NEXTGENDEM/bcs-backend/
elif [ "$(whoami)" == "acurbelo" ] ; then
  cd ~/Proyectos/NEXTGENDEM/bcs-backend/
elif [ "$(whoami)" == "pmoreno" ] ; then
  cd ~/Proyectos/NEXTGENDEM/bcs-backend/
elif [ "$(whoami)" == "dreyes" ] ; then
  cd ~/Proyectos/NEXTGENDEM/bcs-backend/
fi

# CELERY
# Worker AND Beat (only for development; NOT for production -use Supervisor and two separate processes-)
celery -A biobarcoding.tasks.definitions worker --beat --loglevel=info
