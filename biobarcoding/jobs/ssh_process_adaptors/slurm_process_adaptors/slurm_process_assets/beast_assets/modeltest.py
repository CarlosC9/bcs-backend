import os

MODELTEST_FOLDER = "/home/dreyes/beast_pipeline"

def get_substitution_model(filename: str):
    bic_model_first_line = "Best model according to BIC"
    model = None
    invProp = None
    gamma = None
    with open(filename, "r") as f:
        first_line_found = False
        model_found = False
        inv_found = False
        gamma_found = False
        for l in f:
            if l.strip() == bic_model_first_line:
                first_line_found = True
            if first_line_found:
                if l.startswith("Model"):
                    if "TrN" in l.split()[1]:
                        model = "TN93"
                    else:
                        model = "HKY"
                    model_found = True
                elif l.startswith("Inv."):
                    if l.split(":")[1].strip() != "-":
                        invProp = float(l.split(":")[1].strip())
                    inv_found = True
                elif l.startswith("Gamma"):
                    if l.split(":")[1].strip() != "-":
                        gamma = float(l.split(":")[1].strip())
                    gamma_found = True
                if model_found and inv_found and gamma_found:
                    return {"model": model, "invProp": invProp, "gamma": gamma}
    return "error"

def get_substitution_models(alignments_filenames: [], threads: str):
    sub_models = {}
    for aln_f in alignments_filenames:
        if not os.path.exists(f"{aln_f.split('.')[0]}_modeltest.out"):
            os.system(f"{os.getenv('BEAST_DEPENDENCIES_PATH')}/modeltest-ng-mpi -p {threads} --disable-checkpoint -m HKY,TrN -i {aln_f} -o {aln_f.split('.')[0]}_modeltest")
        sub_model = get_substitution_model(f"{aln_f.split('.')[0]}_modeltest.out")
        sub_models[aln_f.split(".")[0]] = sub_model
    return sub_models


if __name__ == '__main__':
    import sys
    sub_model = get_substitution_models(sys.argv[1].split(","), sys.argv[2])
    