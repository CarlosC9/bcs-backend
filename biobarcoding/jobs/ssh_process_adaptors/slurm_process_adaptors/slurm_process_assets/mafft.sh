#!/bin/bash

export n=`expr $ntasks \* $cpus_per_task`
export MAFFT_N_THREADS_PER_PROCESS="$cpus_per_task"
export MAFFT_MPIRUN="mpirun -n $n --map-by ppr:$cpus_per_task:node -bind-to none --oversubscribe"

if [[ $SLURM_CLUSTER_NAME == teidehpc ]]
then
  source /etc/profile.d/profile.modules.sh
  module load gcc/10.2.0
  module load openmpi/3.1.4/gcc
  module load mafft/7.487/openmpi
fi

if [[ $algorithm = "--6merpair" ]]
then
  srun mafft --mpi --large $algorithm --retree $retree $parttree --maxiterate $maxiterate $input_filename > $output_filename
else
  srun mafft --mpi --large $algorithm --maxiterate $maxiterate $input_filename > $output_filename
fi