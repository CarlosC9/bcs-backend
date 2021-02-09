#!/bin/bash
# Assumes the other Docker Compose containers are up
# To ensure they are running:
#  - login into balder
#  - cd /mnt/datos/bcs_projects/sys/bc_ecosys/dev_balder
#  - docker-compose logs -f (the log should appear)
#  - docker-compose ps (a list of services should appear)
#  - if no service appears, maybe all are down. To start: cd .. && ./start.sh dev_balder/environment.sh
#  - if a shutdown is wanted: docker-compose down

# Move "bcs-backend" files to BALDER
rsync -avhz --delete --exclude ".*/" /home/rnebot/GoogleDrive/AA_NEXTGENDEM/bcs-backend/ $(whoami)@10.141.150.137:/mnt/datos/bcs_projects/backend
# Build, stop and start "backend"
ssh $(whoami)@10.141.150.137 'cd /mnt/datos/bcs_projects/sys/bc_ecosys/dev_balder && docker-compose build backend && docker-compose stop backend && docker-compose up -d --no-deps backend'
