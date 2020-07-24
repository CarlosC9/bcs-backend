#!/bin/bash
# Prepare Celery for development of "bcs-backend"
docker run -name redis --rm -d -p 6379:6379 redis
cd ~/GoogleDrive/AA_NEXTGENDEM/bcs-backend/
# Worker AND Beat (only for development; NOT for production -use Supervisor and two separate processes-)
celery -A biobarcoding.tasks.main worker --beat --loglevel=info