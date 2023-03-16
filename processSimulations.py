import pytraj as pt
import numpy as np 
import pandas as pd
import matplotlib.pylab as plt

dataFile = "traj3.csv"
#fileToProcess = "system_reduced_protein.pdb"

fileToProcess = "../trajs7_ad/system_reduced_protein.pdb"
dataFile = "traj7.csv"

# inputs 
protLen=202
# firstAtom = "BB  MET B   1"
mask = ":MET@BB" # i can't get resid 1 to work correctly, so this is a workaround



# find number of structures

d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def ProcessTraj(fileName=""): 
  traj = pt.iterload(fileToProcess)                 
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
    daHisto,binEdges = np.histogram(data, bins=10,range=[10,20],density=True)
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
  feature = "fasta"
  #nEle = len( df.index ) 
  vals = df[feature]
  nRes = int(len(vals[0]))
  scores = [ x.count("E")/nRes for x in vals ]
  df['negative']=scores

import re 
def stringArToAr(stringAr):
  val = re.sub('\[\s*',"",stringAr)
  val = re.sub('\s*\]',"",val)
  val = re.sub('\n',"",val)
  print(val)
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
task="GenerateData"
task="ProcessData"

if task is "GenerateData":
  nStruct,traj = ProcessTraj()           
  nStruct=30 
  print("Processing %d"%nStruct)
  copyData = GetTrajData(traj,nStruct = nStruct)
  df = pd.DataFrame.from_dict(copyData) 

  # should do pickle eventually) 
  df.to_csv(dataFile) 

elif task is "ProcessData":
  inputFile = dataFile             
  dfa = pd.read_csv( inputFile )              

  ScoreFasta(dfa)

  out = dataFile.replace('.csv',"_scored.csv") 
  dfa.to_csv(out)                              





  
