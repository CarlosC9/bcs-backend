#!/bin/bash
# Prepare Celery for development of "bcs-backend"
if [ ! "$(docker ps -q -f name=redis)" ]; then
  docker run -name redis --rm -d -p 6379:6379 redis
fi

if ["$(whoami)" == "rnebot"]; then
  cd ~/GoogleDrive/AA_NEXTGENDEM/bcs-backend/
elif [ "$(whoami)" == "acurbelo" ]; then
  cd .
fi

# Worker AND Beat (only for development; NOT for production -use Supervisor and two separate processes-)
celery -A biobarcoding.tasks.definitions worker --beat --loglevel=info