#!/bin/bash

if [[ $SLURM_CLUSTER_NAME == teidehpc ]]
then
  source /etc/profile.d/profile.modules.sh
  module load gcc/10.2.0
  module load openmpi/3.1.4/gcc
  module load mafft/7.487/openmpi
fi

if [[ $algorithm = "--6merpair" ]]
then
  srun --nodes=$n --ntasks-per-node=$cpus_per_task mafft --mpi --large $algorithm --thread $cpus_per_task --retree $retree $parttree --maxiterate $maxiterate $input_filename > $output_filename
else
  srun --nodes=$n --ntasks-per-node=$cpus_per_task mafft --mpi --large $algorithm --thread $cpus_per_task --maxiterate $maxiterate $input_filename > $output_filename
fi