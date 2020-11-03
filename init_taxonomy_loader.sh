#!/bin/bash

cd $(dirname $0)/biobarcoding/services/perl_scripts/

wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz.md5
while [ "$(md5sum taxdump.tar.gz)" != "$(cat taxdump.tar.gz.md5)" ]
do
  wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
  wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz.md5
done

tar -zxvf taxdump.tar.gz names.dmp nodes.dmp
rm taxdump.tar.gz taxdump.tar.gz.md5

perl ./load_taxonomy_cvterms.pl -H localhost -D postgres -u postgres -d Pg -p postgres
