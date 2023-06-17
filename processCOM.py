import concurrent.futures
import json
import os
import re
import sys
import time
from itertools import repeat

import numpy as np
import pandas as pd
import pytraj as pt

## specific to the protein

protLen = 202
mask = ":MET@BB"  ## i can't get resid 1 to work correctly, so this is a workaround
nStruct = 200

d = {
    "CYS": "C",
    "ASP": "D",
    "SER": "S",
    "GLN": "Q",
    "LYS": "K",
    "ILE": "I",
    "PRO": "P",
    "THR": "T",
    "PHE": "F",
    "ASN": "N",
    "GLY": "G",
    "HIS": "H",
    "LEU": "L",
    "ARG": "R",
    "TRP": "W",
    "ALA": "A",
    "VAL": "V",
    "GLU": "E",
    "TYR": "Y",
    "MET": "M",
}


def LoadTraj(caseToProcess):
    ## load

    if "pdb" not in caseToProcess:
        trajName = caseToProcess + ".xtc"
        protName = caseToProcess + ".pdb"
        traj = pt.iterload(trajName, protName)
    else:
        traj = pt.iterload(caseToProcess)

    ## get structure info

    fr_i = traj[0, mask]
    l = pt.get_coordinates(fr_i)
    nStruct = np.shape(l)[1]
    print("Finding %d structures" % nStruct)

    return nStruct, traj


def GetTrajData(traj, nStruct):
    ## prints out summary of the trajectory from the simulation
    print(traj)

    pt._verbose()
    copies = []

    ## goes through each copy of the trajectory sequentially for a maximum of 200 according to documentation

    start = ((nStruct-1) * protLen) + 1
    fin = start + protLen - 1
    mask = "@%d-%d" % (start, fin)  ############ find out what the d's mean

    ## superpose
    pt.superpose(traj, mask=mask)

    com = pt.center_of_mass(traj, mask=mask)

    return com.tolist()


def doit(nStruct, mode=None, case=None):
    if "trajs3" in case:
        caseToProcess = os.path.join(case, "system_reduced_all.pdb")
        dataSeries = "./traj3COMTest.json"
        # dataSeries = "./traj3COM.csv"
        # dataSeries = "./data/traj3Series.json"

    elif "trajs7" in case:
        caseToProcess = os.path.join(case, "system_reduced_all_MDtraj.pdb")
        dataSeries = "./traj7COMTest.json"
        # dataSeries = "./traj7COM.csv"
        # dataSeries = "./data/traj7Series.json"
        
    else:
        raise RuntimeError("dunno this case")

    ## inputs

    if mode is "generation":
        print("Generating data from trajs")
        nStructPossible, traj = LoadTraj(caseToProcess)
        nStruct = np.min([nStruct, nStructPossible])
        print("Processing %d" % nStruct)

        ## saving series data

        with concurrent.futures.ProcessPoolExecutor(max_workers=25) as executor:
            indices = [i for i in range(1, nStruct+1)]
            results = executor.map(GetTrajData, repeat(traj), indices)

        seriesData = {i: result for i, result in enumerate(results)} 
        # seriesData = [result for result in results]

        # seriesData = np.array(seriesData)

        # seriesData = pd.DataFrame(seriesData)

        # seriesData.to_csv(dataSeries)
        
        with open(dataSeries, 'w') as file:
            json.dump(seriesData, file)


    elif mode is "postprocess":
        print("Postprocessing data from trajs")
        inputFile = dataFile

        val = re.sub("traj", "", case)

        out = dataFile.replace(".csv", "_scored.csv")
        dfa.to_csv(out)
    else:
        raise RuntimeError("mode not understood")


def helpmsg():
    scriptName = sys.argv[0]
    msg = """
Purpose: 
 
Usage:
"""
    msg += "  %s -validation" % (scriptName)
    msg += """

Notes:

"""
    return msg


## main code running script when run from commandline

mode = "generation"
if __name__ == "__main__":
    msg = helpmsg()
    remap = "none"

    if len(sys.argv) < 2:
        raise RuntimeError(msg)

    # Loops over each argument in the command line
    for i, arg in enumerate(sys.argv):
        # print(arg)
        # calls 'doit' with the next argument following the argument '-validation'
        if arg == "-generation":
            mode = "generation"
        if arg == "-postprocess":
            mode = "postprocess"
        if arg == "-nstruct":
            nStruct = int(sys.argv[i + 1])
        if arg == "-case":
            arg1 = sys.argv[i + 1]
            start = time.perf_counter()
            doit(nStruct, mode=mode, case=arg1)
            end = time.perf_counter()
            print(f"Time to run = {end-start}")
            quit()

    raise RuntimeError("Arguments not understood")