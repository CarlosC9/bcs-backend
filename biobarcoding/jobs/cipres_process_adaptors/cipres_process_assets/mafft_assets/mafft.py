import os
import sys

distanceMetricAlgorithm = {"--6merpair": 0, "--globalpair": 1, "--localpair": 2, "--genafpair": 3, "--fastapair": 4}
partTreeMetric = {"--parttree": 0, "--dpparttree": 1}

if __name__ == "__main__":
    infile = os.path.join(os.path.dirname(__file__), sys.argv[8])
    inputParams = {"infile_": infile}
    print(sys.argv)
    vParams = {
        "tool": "MAFFT_XSEDE",
        "runtime_": float(sys.argv[10]) / 60,
        "auto_analysis_": 0,
        "configure_analysis_": 1,
        "distanceMetric_": distanceMetricAlgorithm[sys.argv[6]],
        "iterativeRefinements_": int(sys.argv[7]),
    }

    if sys.argv[6] == "--6merpair":
        vParams["retrees_"] = int(sys.argv[11])
        if len(sys.argv) > 12:
            vParams["usePartTree_"] = 1
            vParams["partTreeMetric_"] = partTreeMetric[sys.argv[12]]
    print("-----------INPUT PARAMS--------------")
    print(inputParams)
    print("-----------V PARAMS--------------")
    print(vParams)
    try:
        from submit_cipres import submit_cipres
        submit_cipres(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], vParams, inputParams)
    except Exception as e:
        print("Before submit error")
        print(e)

