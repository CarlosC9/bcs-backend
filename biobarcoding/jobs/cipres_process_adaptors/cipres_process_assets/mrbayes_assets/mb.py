import os
import sys
import json

argument_names = ['start_phylotree', 'nst', 'rates', 'taxons_select', 'ngen', 'nchains', 'samplefreq', 'filename', 'burninfrac']

def write_template(argv):
    #Read alignment
    with open('aln.nexus', 'r') as f:
        alignment = f.read()

    if argv[0] == "true":  #if guide tree is used
        with open('start_phylotree.nexus', 'r') as f:
            tree = f.read()
            tree = tree.replace("#NEXUS", "")

    for i in range(1, len(argv)):
        template = template.replace(f"${argument_names[i]}", argv[i])

    final_script = f"{alignment}\n{tree}\n{template}"
    print(final_script)
    with open('mb_batch.nex', 'w') as f:
        f.write(final_script)

if __name__ == "__main__":
    write_template(sys.argv[6:])
    infile = os.path.join(os.path.realpath(__file__), "mb_batch.nex")
    inputParams = {"infile_": infile}
    vParams = {
        "toolId": "MRBAYES_XSEDE",
        "runtime_": 0.5,
        "mrbayesblockquery_": 1,
        "nchains_specified_": sys.argv[11],
        "nruns_specified_": 4,
        "set_beagle_params_": 1,
    }
    cmd = f"python3 submit_cipres.py {' '.join(sys.argv[1:6])} \"{json.dumps(vParams)}\" \"{json.dumps(inputParams)}\""
    sys
