# !/bin/bash
# set -o errexit
set -o nounset

# A bash script that gets the taxon IDs from NCBI for a given list of taxon names
# INPUT: A file with the list of taxon names
# OUTPUT: A file with the list of taxon IDs
# WATCH_OUT: OUTPUT contains !Lineage !Show organism modifiers

function init_sh() {
  readonly input_file='biota2ncbi_input';
  readonly tmp_file='curl_response';
  readonly output_file='biota2ncbi_output';
  readonly taxon_finder='wwwtax.cgi?mode=';
  # SHORT_LINEAGE: QUERY with lin=s
  readonly ncbi_query='https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?lin=s&name=';
  echo ' ## VARIABLES DECLARED'
}

function run_sh() {
  echo ' ## RUNNING'
  echo > "$tmp_file"
  while read name; do
    echo " >> $name"
    query=$(echo "$ncbi_query$name" | sed 's/ /%20/;')
    echo "$query"
    curl -s "$query" | grep "$taxon_finder" >> "$tmp_file"
  done < <(cut -d ' ' -f 1-2 "$input_file")
  # FILTERING: sed 's/id=/\n/g' "$tmp_file" | grep -f <taxon_type> | ...
  sed 's/id=/\n/g' "$tmp_file" | sed 's/^\([0-9]*\).*/\1/; /^$/d' | sort -h -u > "$output_file"
  rm "$tmp_file"
  echo ' ## DONE'
}

init_sh
run_sh
exit 0
