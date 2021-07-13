#!/bin/bash

TEMPLATE_FILENAME=paup_parsimony.txt
INPUTS_FILENAME=ngd_paup_parsimony.txt

check_error()
{
  if [ $1 -ne 0 ]; then
    echo "Failure. Exit Status: $1"
    exit $1
  fi
}
echo $#
if [[ $# -eq 15 ]]
then
  outRoot=$1
  gapMode=$2
  addseq=$3
  swap=$4
  hold=$5
  consensus_tree_type=$6
  le50=$7
  percent=$8
  n_replicas=$9
  search=${10}
  method=${11}
  enforce_converse=${12}
  taxset_paup=${13}
  sets=${14}
  assumptions=${15}
  python3 paup_parsimony_params.py $TEMPLATE_FILENAME $INPUTS_FILENAME $outRoot $gapMode $addseq $swap $hold "$consensus_tree_type" $le50 $percent $search $n_replicas $method "$enforce_converse" "$taxset_paup" "$sets" "$assumptions"
  paup ngd_paup_parsimony.txt
else
  echo "Something went wrong with the arguments." >&2
  exit 1
fi

check_error $?