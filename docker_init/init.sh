#!/bin/bash
requires=("perl" "libdbi-perl" "libgo-perl" "libpq-dev" "cpanminus" "bioperl")
dpkg -s "${requires[@]}" >/dev/null 2>&1
if [ ! $? ]; then
  sudo apt -y install perl libdbi-perl libgo-perl libpq-dev cpanminus bioperl
  sudo cpan install Bundle::GMOD
  sudo cpanm GO::Utils DBIx::DBStag DBIx::DBSchema DBD::Pg
  # Chado Perl "Bio" module
  wget https://github.com/GMOD/Chado/archive/master.zip -O /tmp/chado_src.zip
  unzip -o /tmp/chado_src.zip -d /tmp/chado_project
  sudo cp -r /tmp/chado_project/Chado-master/chado/lib/Bio /usr/local/share/perl/5.30.0/
  rm /tmp/chado_src.zip /tmp/chado_project -r
  sudo cpan Bio::Chado::Schema
fi
python3 python_scripts/edam_insertion.py -h localhost
./init_taxonomy_loader.sh
