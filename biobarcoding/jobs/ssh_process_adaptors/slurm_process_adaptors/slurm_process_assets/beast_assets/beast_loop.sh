#!/bin/bash

list_descendants ()
{
  local children=$(ps -o pid= --ppid "$1")

  for pid in $children
  do
    list_descendants "$pid"
  done

  echo "kill ${children}"
}

nohup bash -c "$BEAST_DEPENDENCIES_PATH/beast/bin/beast -beagle_CPU -threads $2 $1" &
pid=$!
pgid=$(ps -p $pid -o pgid=)
echo "---PID--- ${pid}"
echo "---PGID--- ${pgid}"
echo $my_dir/tracerer.R
converged=TRUE
ps --pid $pid > /dev/null
counter=0
while [ $? -eq 0 ] && [ $counter -lt 2 ]; do
	converged=$(Rscript ../../tracerer.R)
	if [ $converged == FALSE ]; then
		counter=$(($counter+1))
	else
		counter=0
	fi
	sleep 10
	ps --pid $pid > /dev/null
done
kill $(list_descendants $$)
