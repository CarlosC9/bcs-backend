import os

PATH_SAMPLING_HEADER = """<run spec='beast.inference.PathSampler' nrOfSteps='100' chainLength='1000000' alpha='0.3' rootdir ='./tmp/step' burnInPercentage='50' preBurnin='1000000' deleteOldLogs='true'>
	cd $(dir)
 	java -cp $(java.class.path) beast.app.beastapp.BeastMain $(resume/overwrite) -java -seed $(seed) beast.xml\n"""
MCMC_HEADER = '<mcmc id="mcmc" spec="MCMC" chainLength="100000000">\n'


def generate_path_sampling(filename: str):
    path_sampling_content = ""
    with open(filename, 'r') as f:
        for l in f:
            if l.strip().startswith("<run"):
                path_sampling_content += PATH_SAMPLING_HEADER
                path_sampling_content += MCMC_HEADER
            elif  "</run>" in l:
                path_sampling_content += l.replace("</run>", "\n</mcmc>\n</run>")
            else:
                path_sampling_content += l
    
    if path_sampling_content:
        with open(os.path.basename(filename).split(".")[0] + "_path_sampling.xml", "w") as f:
            f.write(path_sampling_content)

if __name__ == '__main__':
    import sys
    generate_path_sampling(sys.argv[1])
