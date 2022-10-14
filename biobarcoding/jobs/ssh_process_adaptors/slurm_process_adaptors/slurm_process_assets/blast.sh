#!/bin/bash

if [[ $task == *"megablast"* ]] && [[ $db == "jardin_botanico" ]]; then
  srun blastn -query $query_filename -task $task -db "/home/dreyes/blast/jardin_botanico" -num_threads $cpus_per_task -use_index true -out $output_filename
elif [[ $db == "jardin_botanico" ]]; then
  srun blastn -query $query_filename -task $task -db "/home/dreyes/blast/jardin_botanico" -num_threads $cpus_per_task -out $output_filename
elif [[ $task == *"megablast"* ]] && [[ $db == "nt" ]]; then
  srun blastn -query $query_filename -task $task -db "/mnt/datos/blast_databases/nt/nt" -num_threads $cpus_per_task -use_index true -out $output_filename
elif [[ $db == "nt" ]]; then
  srun blastn -query $query_filename -task $task -db "/mnt/datos/blast_databases/nt/nt" -num_threads $cpus_per_task -out $output_filename
fi