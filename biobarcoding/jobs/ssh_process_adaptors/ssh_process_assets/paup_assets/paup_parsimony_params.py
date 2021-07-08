from sys import argv


def replace_inputs_dict(file_content: str, inputs_dict: dict):
    for key, value in inputs_dict.items():
        file_content = file_content.replace(key, value)
    return file_content


def delete_method_comment(file_content: str, method: str, search: str):
    if search == "fastStep" or search == "BandB":
        file_content = file_content.splitlines()
        new_file_content = []
        for line in file_content:
            if method in line and '/' in line:
                line = line.split('/')[0] + ";"
                line = line.replace(f"[${method} ", "")
            new_file_content.append(line)
        file_content = "\n".join(new_file_content)

    file_content = file_content.replace(f"[${method} ", "")
    file_content = file_content.replace(f" ${method}]", "")
    return file_content


def create_sets_file(sets: str, assumptions: str):
    with open("sets_and_assumptions.txt", "w") as f:
        f.write("#NEXUS\n\n")
        f.write('begin sets' + ';\n')
        f.write(sets.replace("\\\\", "\\"))
        f.write('end;\n\n')
        f.write('begin assumptions' + ';\n')
        f.write(assumptions.replace("\\\\", "\\"))
        f.write('end;')


if __name__ == "__main__":
    template_filename, output_filename, out_root, gap_mode, addseq, swap, hold, consensus_tree_type, le50, percent, n_replicas, search, method, sets, assumptions = argv[1:]
    inputs_dict = {
        "$outRoot": out_root,
        "$gapMode": gap_mode,
        "$addseq": addseq,
        "$swap": swap,
        "$hold": hold,
        "$consensus_tree_type": consensus_tree_type,
        "$le50": f"le50={le50}" if le50 != 'None' else "",
        "$percent": f"percent={percent}" if percent != 'None' else "",
        "$nReplicas": n_replicas,
        "$search": search
    }

    # Generate Paup script with user's parameters
    with open(template_filename, "r") as f:
        file_content = f.read()

    file_content = replace_inputs_dict(file_content, inputs_dict)
    file_content = delete_method_comment(file_content, method, search)

    with open(output_filename, "w") as f:
        f.write(file_content)

    #Create Paup sets and assumptions in sets_and_assumptions.txt file
    create_sets_file(sets, assumptions)
