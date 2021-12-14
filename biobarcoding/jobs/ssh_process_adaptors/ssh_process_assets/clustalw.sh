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


if [[ $# -eq 7 ]]
then
  if [[ $6 -eq 'none' ]]
  then
    clustalw -INFILE="$1" -OUTFILE="$2" -OUTORDER=$3 -TYPE=$4 -OUTPUT=$5
  else
    clustalw -INFILE="$1" -OUTFILE="$2" -OUTORDER=$3 -TYPE=$4 -OUTPUT=$5 -RANGE=$6,$7
  fi
else
  echo "Something went wrong with the arguments." >&2
  exit 1
fi

check_error $?
mv input_dataset.dnd input_dataset.newick