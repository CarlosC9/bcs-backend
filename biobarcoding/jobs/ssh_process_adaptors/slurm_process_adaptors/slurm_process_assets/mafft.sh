#!/bin/bash

export n=`expr $ntasks \* $cpus_per_task`
export MAFFT_N_THREADS_PER_PROCESS="$cpus_per_task"
export MAFFT_MPIRUN="mpirun -n $n -N $cpus_per_task -bind-to none --oversubscribe" #ntasks processes with 1 thread

if [[ $algorithm = "--6merpair" ]]
then
  srun mafft --mpi --large $algorithm --retree $retree $parttree --maxiterate $maxiterate $input_filename > $output_filename
else
  srun mafft --mpi --large $algorithm --maxiterate $maxiterate $input_filename > $output_filename
fi