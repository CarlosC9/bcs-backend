from sys import argv
import dendropy

if __name__ == "__main__":
    nexus_file = argv[1]
    tree = dendropy.Tree.get(path=nexus_file, schema="nexus")
    tree.write(path=nexus_file.rsplit(".", 1)[0] + ".nexus", schema="nexus", suppress_taxa_blocks=True)
