#!/bin/bash

if [[ $algorithm = "--6merpair" ]]
then
  python3 mafft.py $base_url $username $password $appID $app_name $algorithm $maxiterate $input_filename $cpus_per_task $time $retree $parttree
else
  python3 mafft.py $base_url $username $password $appID $app_name $algorithm $maxiterate $input_filename $cpus_per_task $time
fi