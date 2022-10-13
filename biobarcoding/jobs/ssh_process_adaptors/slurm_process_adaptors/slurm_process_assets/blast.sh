#!/bin/bash

if [[ $task == *"megablast"* ]]; then
  srun blastn -query $query_filename -task $task -db "/home/dreyes/blast/$db" -num_threads $cpus_per_task -use_index true -out $output_filename
else
  srun blastn -query $query_filename -task $task -db "/mnt/datos/blast_databases/$db" -num_threads $cpus_per_task -out $output_filename
fi