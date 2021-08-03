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

chmod +x clustalw.sh
chmod +x paup_parsimony.sh
chmod +x pda.sh

if [[ $# -eq 22 ]]
then
  ./clustalw.sh "$1" "$2" $3 $4 $5
  check_error $?
  python3 fasta2nexus.py clustal.fa alignment.txt . -
  check_error $?
  ./paup_parsimony.sh $6 $7 $8 $9 ${10} "${11}" ${12} ${13} ${14} ${15} ${16} "${17}" "${18}" "${19}" "${20}"
  check_error $?
  cp *_consensus.tre consensus.tre
  check_error $?
  python3 nexus2newick.py consensus.tre consensus.newick
  check_error $?
  ./pda.sh ${21} ${22}
  check_error $?
elif [[ $# -eq 24 ]]
then
  ./clustalw.sh "$1" "$2" $3 $4 $5 $6 $7
  check_error $?
  python3 fasta2nexus.py clustal.fa alignment.txt . -
  check_error $?
  ./paup_parsimony.sh $8 $9 ${10} ${11} ${12} "${13}" ${14} ${15} ${16} ${17} ${18} "${19}" "${20}" "${21}" "${22}"
  check_error $?
  cp *_consensus.tre consensus.tre
  check_error $?
  python3 nexus2newick.py consensus.tre consensus.newick
  check_error $?
  ./pda.sh ${21} ${22}
  check_error $?
else
  echo "Something went wrong with the arguments." >&2
  exit 1
fi