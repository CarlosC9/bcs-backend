#!/bin/bash
# REDIS
if [ ! "$(docker ps -q -f name=redis)" ]; then
  docker run -name redis --rm -d -p 6379:6379 redis
fi

# PostgreSQL
if [ ! "$(docker ps -q -f name=postgres_devel)" ]; then
  if ["$(whoami)" == "rnebot"]; then
    docker run -n postgres_devel -p 5432:5432 --rm -e POSTGRES_PASSWORD=postgres -v /home/rnebot/DATOS/pg_devel:/var/lib/postgresql/data postgres
  elif ["$(whoami)" == "acurbelo"]; then
    echo docker run command
  fi
fi

# CD
if ["$(whoami)" == "rnebot"]; then
  cd ~/GoogleDrive/AA_NEXTGENDEM/bcs-backend/
elif [ "$(whoami)" == "acurbelo" ]; then
  cd .
fi

# CELERY
# Worker AND Beat (only for development; NOT for production -use Supervisor and two separate processes-)
celery -A biobarcoding.tasks.definitions worker --beat --loglevel=info