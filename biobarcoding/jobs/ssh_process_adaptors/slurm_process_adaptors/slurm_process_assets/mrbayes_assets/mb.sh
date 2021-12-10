#!/bin/bash

python3 mb_template_writer.py  "$input_tree_cmd" $nst $rates $taxons_select $ngen $nchains $samplefreq $filename $burninfrac
srun mpirun -np $ntasks mb mb_batch.nex