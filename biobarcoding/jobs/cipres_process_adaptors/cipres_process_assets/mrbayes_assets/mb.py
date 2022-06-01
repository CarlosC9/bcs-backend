import os
import sys

argument_names = ['nst', 'rates', 'taxons_select', 'ngen', 'nchains', 'samplefreq', 'filename', 'burninfrac']

def write_template(argv):
    #Read alignment
    with open('aln.nexus', 'r') as f:
        alignment = f.read()

    with open('mb_batch_template.nex', 'r') as f:
        template = f.read()

    for i in range(len(argv)):
        template = template.replace(f"${argument_names[i]}", argv[i])

    final_script = f"{alignment}\n{template}"
    print(final_script)
    with open('mb_batch.nex', 'w') as f:
        f.write(final_script)

if __name__ == "__main__":
    write_template(sys.argv[6:14])
    infile = os.path.join(os.path.dirname(__file__), "mb_batch.nex")
    inputParams = {"infile_": infile}
    print(sys.argv)
    vParams = {
        "tool": "MRBAYES_XSEDE",
        "runtime_": float(sys.argv[15]) / 60,
        "mrbayesblockquery_": 1,
        "nchains_specified_": sys.argv[10],
        "nruns_specified_": 1,
        "set_beagle_params_": 1,
    }
    print("-----------INPUT PARAMS--------------")
    print(inputParams)
    print("-----------V PARAMS--------------")
    print(vParams)
    from submit_cipres import submit_cipres
    submit_cipres(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], vParams, inputParams)
