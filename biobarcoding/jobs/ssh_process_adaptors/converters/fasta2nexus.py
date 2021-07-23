from bioconvert.fasta2nexus import FASTA2NEXUS
from sys import argv

if __name__ == "__main__":
    fasta_file, nexus_file, missing, gap = argv[1:]
    FASTA2NEXUS(fasta_file, nexus_file)()
    with open(nexus_file, "r") as f:
        lines = f.readlines()
        print(lines)
        lines[3] = f"format datatype=dna missing={missing} gap={gap};\n"
    with open(nexus_file, "w") as f:
        f.writelines(lines)