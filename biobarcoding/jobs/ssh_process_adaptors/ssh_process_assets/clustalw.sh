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

if [[ $# -eq 5 ]]
then
  #./clustalw.sh "matK_25taxones_Netgendem_SINalinear.fasta" "clustalw.aln" ALIGNED DNA
  clustalw -INFILE="$1" -OUTFILE="$2" -OUTORDER=$3 -TYPE=$4 -OUTPUT=$5
elif [[ $# -eq 7 ]]
then
  clustalw -INFILE="$1" -OUTFILE="$2" -OUTORDER=$3 -TYPE=$4 -OUTPUT=$5 -RANGE=$6,$7
else
  echo "Something went wrong with the arguments." >&2
  exit 1
fi

check_error $?