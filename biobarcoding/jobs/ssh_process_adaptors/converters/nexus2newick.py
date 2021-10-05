from sys import argv

from bioconvert.nexus2newick import NEXUS2NEWICK

if __name__ == "__main__":
    fasta_file, nexus_file = argv[1:]
    NEXUS2NEWICK(fasta_file, nexus_file)()
