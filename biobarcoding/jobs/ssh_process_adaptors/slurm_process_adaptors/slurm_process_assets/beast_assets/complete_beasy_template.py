import os
import re
from fasta2nexus import convert_fasta2nexus
import modeltest


def create_imports_str(alignments_dict: []):
    import_template = f"import {os.getenv('my_dir')}/$alignment;\n"
    import_str = ""
    for a in alignments_dict.keys():
        if not os.path.exists(f"{a}.nex"):
            convert_fasta2nexus(f"{a}.fasta", f"{a}.nex", "?", "-")  # convert fasta
        import_str += import_template.replace("$alignment", f"{a}.nex")
    return import_str


def create_alignments_str(alignments_dict: []):
    alignment_template = "{$alignment}"
    return alignment_template.replace("$alignment", ", ".join(list(alignments_dict.keys())))


def create_sub_model_str(region_id, sub_model):
    model_str = f"\nuse substModel {{{region_id}}} = {sub_model['model']};\n"
    if sub_model['gamma']:
        model_str += f"""set siteModel@gammaCategoryCount {{{region_id}}} = 5;\nset shape {{{region_id}}} = {sub_model['gamma']}; 
        set shape@estimate {{{region_id}}} = true;\n"""
    else:
        model_str += f"set siteModel@gammaCategoryCount {{{region_id}}} = 0;\n"
    if sub_model['invProp']:
        model_str += f"""set proportionInvariant {{{region_id}}} = {sub_model['invProp']};
        set proportionInvariant@estimate {{{region_id}}} = true;\n"""
    return model_str


def create_substitution_models_str(alignments_dict: [], threads):
    alignments_filenames = []
    for a in alignments_dict.keys():
        alignments_filenames.append(f"{a}.fasta")
    sub_models = modeltest.get_substitution_models(alignments_filenames, threads)
    models_strings = ""
    for region_id, sub_model in sub_models.items():
        models_strings += create_sub_model_str(region_id, sub_model)
    return models_strings


def create_clock_models_str(alignments_dict: []):
    clock_models_template = """use branchRateModel {$region_id} = RelaxedClockLogNormal;
    set clock.rate@estimate {$region_id} = true;\n"""
    clock_models_str = ""
    for a in alignments_dict.keys():
        clock_models_str += clock_models_template.replace("$region_id", a)
    return clock_models_str


def create_monophyly_str(monophyly_taxons: {}):
    monophyly_str = ""
    if monophyly_taxons:
        for key, taxons in monophyly_taxons.items():
            monophyly_str += f"""taxonset {key} = {taxons};
            add MRCAPrior({key}, Uniform(lower=0,upper=100000000), treepartition);
            set monophyletic[{key}.prior] = true;\n"""
    return monophyly_str


def dict_input2dict(dict_input: str):

    new_dict = {}
    if dict_input != "":
        result = re.findall(r"[^{}:,]+", dict_input)
        for i in range(0, len(result), 2):
            new_dict[result[i].strip()] = result[i+1].strip()
    return new_dict


if __name__ == '__main__':
    import sys
    alignments_dict = dict_input2dict(sys.argv[1])
    method_str = sys.argv[2]
    threads = sys.argv[3]
    monophyly_taxons = dict_input2dict(sys.argv[4])
    imports_str = create_imports_str(alignments_dict)
    alignments_str = create_alignments_str(alignments_dict)
    first_alignment_str = list(alignments_dict.keys())[0]
    substitution_models_str = create_substitution_models_str(alignments_dict, threads)
    clock_models_str = create_clock_models_str(alignments_dict)
    monophyly_str = create_monophyly_str(monophyly_taxons)
    with open("beast_template.bea", "r") as f:
        template = f.read()
    if len(alignments_dict.keys()) > 1:
        template = template.replace("$linktree", f"link tree {alignments_str};")
    else:
        template = template.replace("$linktree", "")

    template = template.replace("$template",
                                f"template {os.getenv('BEAST_DEPENDENCIES_PATH')}/beast/templates/myTemplate.xml")
    template = template.replace("$method", method_str)
    template = template.replace("$imports", imports_str)
    template = template.replace("$alignments", alignments_str)
    template = template.replace("$first_alignment", first_alignment_str)
    template = template.replace("$substitution_models", substitution_models_str)
    template = template.replace("$clock_models", clock_models_str)
    template = template.replace("$monophyly", monophyly_str)
    template = template.replace("$threads", threads)
    with open(f"beast_{method_str}.bea", "w") as f:
        f.write(template)
