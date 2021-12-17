#!/bin/bash

python3 mb_template_writer.py  "$input_tree_cmd" $nst $rates $taxons_select $ngen $nchains $samplefreq $filename $burninfrac
srun mpirun -np $ntasks mb mb_batch.nex
for f in *.t; do
  mv -- "$f" "${f%.t}_phylotrees.nexus"
done
for f in *.p; do
  mv -- "$f" "${f%.p}_substitution_model_parameters.tsv"
done
for f in *.mcmc; do
  mv -- "$f" "${f%.mcmc}_convergence_diagnostics.tsv"
done