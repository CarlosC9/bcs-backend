import os
import sys
dirname = os.path.dirname(os.path.abspath(__file__))
print("Adding line")

with open(os.path.join(dirname, "myfile.txt"), "a") as f:
    f.write(f'\n{sys.argv[1]}')