from sys import argv

from bioconvert.fasta2nexus import FASTA2NEXUS


def convert_fasta2nexus(fasta_file, nexus_file, missing, gap):
    FASTA2NEXUS(fasta_file, nexus_file)()
    with open(nexus_file, "r") as f:
        lines = f.readlines()
        lines[3] = f"format datatype=dna missing={missing} gap={gap};\n"
    with open(nexus_file, "w") as f:
        f.writelines(lines)


if __name__ == "__main__":
    fasta_file, nexus_file, missing, gap = argv[1:]
    convert_fasta2nexus(fasta_file, nexus_file, missing, gap)
