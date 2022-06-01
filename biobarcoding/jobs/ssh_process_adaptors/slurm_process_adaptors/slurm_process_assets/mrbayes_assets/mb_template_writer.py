
import sys

argument_names = ['nst', 'rates', 'taxons_select', 'ngen', 'nchains', 'samplefreq', 'filename', 'burninfrac']

def main(argv):
    with open('mb_batch_template.nex', 'r') as f:
        template = f.read()
    
    for i in range(len(argv)):
        template = template.replace(f"${argument_names[i]}", argv[i])
    print(template)
    with open('mb_batch.nex', 'w') as f:
        f.write(template)

if __name__ == "__main__":
    main(sys.argv[1:])