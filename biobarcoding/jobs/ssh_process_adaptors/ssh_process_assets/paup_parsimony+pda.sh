#!/bin/bash
# The #!/bin/bash is necessary to use built-in bash commands with ssh (ssh is not bash)

#Checks if the ouptut of a command is an error. If it is then it exits
check_error()
{
  if [ $1 -ne 0 ]; then
    echo "Failure. Exit Status: $1"
    exit $1
  fi
}

chmod +x paup_parsimony.sh
chmod +x pda.sh
export PATH=/usr/local/go/bin:$PATH
if [[ $# -eq 17 ]]
then
  ./paup_parsimony.sh $1 $2 $3 $4 $5 "$6" $7 $8 $9 ${10} ${11} "${12}" "${13}" "${14}" "${15}"
  check_error $?
  cp *_consensus.nexus consensus.nexus
  check_error $?
  python3 nexus2newick.py consensus.nexus phylotree.newick
  check_error $?
  ./pda.sh ${16} ${17}
  check_error $?
else
  echo "Something went wrong with the arguments." >&2
  exit 1
fi




