import pytraj as pt
import numpy as np 
import pandas as pd
import matplotlib.pylab as plt

##
## Inputs 
##

protLen=202
# firstAtom = "BB  MET B   1"
mask = ":MET@BB" # i can't get resid 1 to work correctly, so this is a workaround
nStruct=200


##
## Functions 
##

# find number of structures
d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def LoadTraj(caseToProcess):
  # load 
  if "pdb" not in caseToProcess:
    trajName = caseToProcess+".xtc"
    protName = caseToProcess+".pdb"
    # xtc file needs associated pdb 
    traj = pt.iterload(trajName,protName)             
  else: 
    traj = pt.iterload(caseToProcess)                 

  # get structure info 
  fr_i = traj[0,mask]
  l = pt.get_coordinates(fr_i); 
  nStruct = np.shape(l)[1]
  print("Finding %d structures"%nStruct)

  return nStruct,traj


def GetFasta(traj,protLen):
    top = traj.top

    # store all residue/indices
    daList=[]
    for residue in top.residues:
      #daList.append(residue)
      daList.append("%d.%s"%(residue.index,residue.name))

    # keep unique (this should 1. keep just one 'tag' for all atoms in a residue and 2. only keep the first protein copy, since all others are identical)
    unique = set( daList ) 

    #print(sorted(unique) )
    l = [ x.split(".")[1] for x in sorted(unique) ]
    fasta = [ d[x] for x in l ]
    fasta = "".join(fasta)

    return fasta





def GetTrajData(traj,nStruct = 2):
  copies = []             
  
  # Assumes all instances are the same 
  # get residue id's 
  fasta = GetFasta(traj,protLen)



  # get per-structure rmsf
  for i in range(nStruct): 
    print(i)
    # define range 
    start=(i*protLen)+1; fin = start+protLen-1
    #mask = "@%d-%d"%(start,fin)
    #fr_i = traj[0,mask]
    #print("%d"%(fin+1)) 
  
    # get atoms for structure
    mask = "@%d-%d"%(start,fin)
    #fr_i = traj[0,mask]

    # superpose
    pt.superpose(traj,mask=mask)
    rmsf = pt.rmsf(traj,mask=mask) 
    # values only, not residues
    np.shape(rmsf) 
    rmsf = rmsf[:,1]
  
    # pt.surf
    # atomiccorr (end to end distance/correlation? ) 
    # ytraj.dihedral_??? miught not work 
    # density? 
    # pca
    # pytraj.watershell
  
    ## compute radgyr
    data = pt.radgyr( traj, mask=mask)
    #plt.plot(data)
    #daHisto,binEdges = np.histogram(data, bins=10,range=[10,20],density=True)
    daHisto,binEdges = np.histogram(data, bins=10,density=True)
    #plt.plot(binEdges,daHisto,label=i)
    #print(binEdges) 
    #plt.plot(binEdges[0:10],daHisto,label=i)

    daHisto=ScoreRg(daHisto,binEdges)
    rmsf =ScoreRMSF(rmsf)         

    container = dict()
    container['copy'] = i      
    container['fasta'] = fasta
    container['binEdges'] = binEdges 
    container['RgHist']=daHisto # if you go back to storing arrays, need to use pickle
    container['RgStart']=data[0]
    container['RgEnd']=data[-1]
    container['RMSF']=rmsf
    #container['salt'] = nearest salt molecules 
    copies.append(container) 

    #plt.plot(rmsf,label=i)
  
  #plt.legend(loc=0)
  #plt.gcf().savefig("test.png") 
  return copies               

def ScoreFasta(df):
  """
  Computes the number of negatively charged a.a. (irrespective of protonation) 
  """
  feature = "fasta"
  #nEle = len( df.index ) 
  vals = df[feature]
  nRes = int(len(vals[0]))

  def tally(vals,aa="E"):
    scores = [ x.count(aa) for x in vals ]
    scores = np.sum(scores)/nRes
    return scores
  nscores = tally(vals,aa="D")
  nscores+= tally(vals,aa="E")
  df['negativeFasta']=nscores

  pscores = tally(vals,aa="K")
  pscores+= tally(vals,aa="R")
  df['positiveFasta']=pscores

def ScoreProtonation(df,pH=None):
  '''
  Protonation state as determined by protonation.ipynb
  HARD CODED
  '''
  pH = int(pH)
  if pH == 3:
      rhoN = -0.052513823529411766 
      rhoP = 0.6387256176470587
  elif pH==7:
      rhoN = -0.621903448275862 
      rhoP = 6.89655172413793e-08
  else:
      raise RuntimeError("pH not understood")

  df['negativepH']=rhoN         
  df['positivepH']=rhoP         
  



import re 
def stringArToAr(stringAr):
  val = re.sub('\[\s*',"",stringAr)
  val = re.sub('\s*\]',"",val)
  val = re.sub('\n',"",val)
  #print(val)
  x = val.split()         
  x = np.asarray(x,dtype='float')
  return x 

def ScoreRg(histo,bins):
  idx = np.argmax(histo) 
  probRg = bins[idx]
  return probRg

def ScoreRMSF(rmsf):
  return np.mean(rmsf)


# get all data 
#nStruct = 10 # 

def doit(mode=None,case=None,nStruct=2):
  #print(case,mode)
  if "traj3" in case:
    caseToProcess = "../trajs3/system_reduced_protein.pdb"
    dataFile = "traj3.csv"
  elif 'traj7' in case:
    caseToProcess = "../trajs7_ad/system_reduced_protein"
    dataFile = "traj7.csv"
  else:
      raise RuntimeError("dunno this case") 

# inputs 

  if mode is "generation":
    print("Generating data from trajs") 
    nStructPossible,traj = LoadTraj(caseToProcess)           
    nStruct = np.min([nStruct,nStructPossible])
    print("Processing %d"%nStruct)
    copyData = GetTrajData(traj,nStruct = nStruct)
    df = pd.DataFrame.from_dict(copyData) 
  
    # should do pickle eventually) 
    print("Printing to ",dataFile) 
    df.to_csv(dataFile) 
  
  elif mode is "postprocess":
    print("Postprocessing data from trajs") 
    inputFile = dataFile             
    dfa = pd.read_csv( inputFile )              
  
    ScoreFasta(dfa)
    val = re.sub('traj',"",case)            
    ScoreProtonation(dfa,pH=val)
  
    out = dataFile.replace('.csv',"_scored.csv") 
    dfa.to_csv(out)                              
  else:
    raise RuntimeError("mode not understood") 
  




  
#!/usr/bin/env python
import sys
##################################
#
# Revisions
#       10.08.10 inception
#
##################################

#
# Message printed when program run without arguments 
#
def helpmsg():
  scriptName= sys.argv[0]
  msg="""
Purpose: 
 
Usage:
"""
  msg+="  %s -validation" % (scriptName)
  msg+="""
  
 
Notes:

"""
  return msg

#
# MAIN routine executed when launching this script from command line 
#
mode="generation"
if __name__ == "__main__":
  import sys
  msg = helpmsg()
  remap = "none"

  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  # Loops over each argument in the command line 
  for i,arg in enumerate(sys.argv):
    #print(arg) 
    # calls 'doit' with the next argument following the argument '-validation'
    if(arg=="-generation"):
      mode="generation"
    if(arg=="-postprocess"):
      mode="postprocess"
    if(arg=="-nstruct"):             
      nStruct=int(sys.argv[i+1]) 
    if(arg=="-case"):
      arg1=sys.argv[i+1] 
      doit(mode=mode,case=arg1,nStruct=nStruct)
      quit()
  





  raise RuntimeError("Arguments not understood")




