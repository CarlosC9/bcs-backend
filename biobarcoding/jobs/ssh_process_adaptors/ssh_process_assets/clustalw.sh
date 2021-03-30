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

if [[ $# -eq 4 ]]
then
  #./clustalw.sh "matK_25taxones_Netgendem_SINalinear.fasta" "clustalw.aln" ALIGNED DNA
  clustalw -INFILE="$1" -OUTFILE="$2" -OUTORDER=$3 -TYPE=$4 -OUTPUT=CLUSTAL
elif [[ $# -eq 6 ]]
then
  #./clustalw.sh "matK_25taxones_Netgendem_SINalinear.fasta" "clustalw.aln" ALIGNED DNA 1 20
  clustalw -INFILE="$1" -OUTFILE="$2" -OUTORDER=$3 -TYPE=$4 -OUTPUT=CLUSTAL -RANGE=$5,$6
else
  echo "Something went wrong with the arguments."
  exit 1
fi

check_error $?