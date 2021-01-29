import os

dirname = os.path.dirname(os.path.abspath(__file__))
print("counting lines")
with open(os.path.join(dirname, "myfile.txt"), "r") as f:
    count = len(f.readlines())
    print(count)
