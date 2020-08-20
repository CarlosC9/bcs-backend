# !/bin/bash
# set -o errexit
set -o nounset

# A bash script that gets the NCBI taxon IDs to populate the phylotree of a given list of taxon IDs
# INPUT: A file with the list of taxon IDs
# OUTPUT: A file with the list of all taxon IDs that make a phylotree

function init_sh() {
  readonly nodes_file='nodes.dmp';
  readonly input_file='orphan_input';
  readonly tmp='tmp';
  readonly parents_file='tmp_parents';
  readonly output_file='orphan_output';
  echo > "$tmp"
  echo > "$parents_file"
  echo > "$output_file"
  echo ' ## VARIABLES DECLARED'
}

function run_sh() {
  echo ' ## RUNNING'
  get_parents "$input_file"
  while [[ 1 -lt $(cat "$parents_file" | wc -l) ]]; do
    get_parents "$parents_file"
    echo "Orphans left: "; cat "$parents_file" | wc -l;
  done
  cat "$parents_file" >> "$output_file"
  sort -h -u "$output_file" > "$tmp"
  mv "$tmp" "$output_file"
  rm "$parents_file"
  echo ' ## DONE'
}

function get_parents() {
  echo ' ## RUNNING: get_parents'
  cat "$1" >> "$output_file"
  sed 's/^\([0-9]*\)$/^\1\\W/' "$1" > "$tmp"
  grep -f "$tmp" "$nodes_file" | cut -d'|' -f2 | tr -d " \t" > "$parents_file"
}

init_sh
run_sh
exit 0
