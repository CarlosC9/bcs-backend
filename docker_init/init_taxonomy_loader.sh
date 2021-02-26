#!/bin/bash

cd $(dirname $0)/../biobarcoding/services/perl_scripts/

if [[ ! -f names.dmp || ! -f nodes.dmp ]]; then
  wget -c ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
  wget -c ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz.md5
  while [ "$(md5sum taxdump.tar.gz)" != "$(cat taxdump.tar.gz.md5)" ]
  do
    echo DIFFERENT
    wget -c ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
    wget -c ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz.md5
  done
  tar -zxvf taxdump.tar.gz names.dmp nodes.dmp
#  rm -f taxdump.tar.gz taxdump.tar.gz.md5
fi

perl ./load_taxonomy_cvterms.pl -H localhost -D postgres -u postgres -d Pg -p postgres
