#Checks if the ouptut of a command is an error. If it is then it exits
check_error()
{
  if [ $1 -ne 0 ]; then
    echo "Failure. Exit Status: $1"
    exit $1
  fi
}

: '
get_output_filename()
{
  # Split str by point. IFS is a special shell variable
  # save original IFS value so we can restore it later
  oIFS="$IFS"
  IFS="."
  declare -a fields=($1)
  IFS="$oIFS"
  unset oIFS
  output_filename="${fields[0]}.$2"
}

get_output_filename "$1" "aln"
'

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
fi

check_error $?