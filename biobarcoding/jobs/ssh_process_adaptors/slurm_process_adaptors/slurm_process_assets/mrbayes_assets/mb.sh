#!/bin/bash

if [[ $SLURM_CLUSTER_NAME == teidehpc ]]
then
  source /etc/profile.d/profile.modules.sh
  module load gcc/10.2.0
  module load openmpi/3.1.4/gcc
  module load beagle/3.1.2/gcc/openmpi
  module load mrbayes/3.2.7/gcc/openmpi
  module load openssl/1.1.1k/gcc
  module load python/3.8.11/gcc
  pip3 install dendropy
fi


python3 mb_template_writer.py $nst $rates $taxons_select $ngen $nchains $samplefreq $filename $burninfrac
srun mpirun --oversubscribe -np $cpus_per_task mb mb_batch.nex
for f in *.t; do
  mv -- "$f" "${f%.t}_phylotrees.t"
done
for f in *.p; do
  mv -- "$f" "${f%.p}_substitution_model_parameters.tsv"
done
for f in *.mcmc; do
  mv -- "$f" "${f%.mcmc}_convergence_diagnostics.tsv"
done

python3 nexus_translate_to_nexus.py "${filename}_phylotrees.t"