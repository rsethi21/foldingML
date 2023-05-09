# foldingML
Original trajectory data are available in 
 /home/pkekeneshuskey/data/molecular_dynamics/ph_binsun
 (from Bin Sun)

## Read this notebook first 
prelim.ipynb - very basic notebook for loading two simulation cases (trajs3 and trajs7) that correspond to an ensemble of proteins simulated at different pHs. Once the files are loaded as pandas dataframes, they are run through a decision-tree classification algorithm. Further details are provided inside the ipynb document. 

trajs3_scored.csv - feature file with classification (isFolded) for proteins simulated at pH 3 
trajs3_scored.csv - feature file with classification (isFolded) for proteins simulated at pH 7 

## To prepare data 
processSimulations.py - this file generates the csv from simulation data. PKH needs to run this, since I have the input simulation files

/home/pkekeneshuskey/data/molecular_dynamics/ph_binsun/foldingML
1. process
python3 processSimulations.py -generate -nstruct 99 -case traj3
python3 processSimulations.py -generate -nstruct 99 -case traj7    
2. postprocess 
python3 processSimulations.py -postprocess -nstruct 99 -case traj3
python3 processSimulations.py -postprocess -nstruct 99 -case traj7    

