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

if [[ $# -eq 2 ]]
then
  pda -root -ts area.nexus $1 $2 phylotree.newick output.pda
  check_error $?
else
  echo "Something went wrong with the arguments." >&2
  exit 1
fi