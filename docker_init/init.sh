#!/bin/bash

sudo apt -y install perl libdbi-perl libgo-perl libpq-dev cpanminus
cpanm GO::Utils DBIx::DBStag DBIx::DBSchema DBD::Pg
python3 python_scripts/edam_insertion.py -h localhost
./init_taxonomy_loader.sh
