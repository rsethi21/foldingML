import argparse
import pandas as pd
import json
import os
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
import numpy as np
import random
from itertools import repeat

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='datafolder with csvs', required=True)

def access(dataFolder):
    
    indexList = []
    copiesList = []

    for item in os.listdir(dataFolder):
        if ".json" in item:
            index = item[item.find("Series")-1]
            with open(os.path.join(dataFolder, item), "r") as file:
                data = json.load(file)
                copies = [pd.DataFrame.from_dict(dictionary) for dictionary in data]
                copiesList = copiesList + copies
                indexList = indexList + [f"{index}-copy{_+1}" for _ in range(len(copies))]
  
    minimum = min([len(copy) for copy in copiesList])
    
    for i, copy in enumerate(copiesList):
        copiesList[i] = copy.iloc[0:minimum,:]
        copiesList[i]['marker'] = [indexList[i] for _ in range(minimum)]

    random.shuffle(copiesList)

    return copiesList

def combine(lst):
    
    return pd.concat(lst)

def normalize(dataframe):
   
    dfToNorm = dataframe.loc[:,dataframe.columns != 'marker']
    scaler = MinMaxScaler()
    scaler.fit(dfToNorm)
    data = scaler.transform(dfToNorm)

    dataframe[dfToNorm.columns] = data

    return dataframe, scaler

def format(dataframe, lag = 30, delay = 0, future = 5):

    numSamples = len(dataframe)
    numSplits = numSamples - lag - delay - future

    input_indices = np.array([range(i, i+lag+1) for i in range(numSplits)])
    output_indices = np.array([range(i+lag+delay, i+lag+delay+future+1) for i in range(numSplits)])

    return input_indices, output_indices

def save():
    return

if __name__ == '__main__':

    args = parser.parse_args()

    dfs = access(args.input)

    dfMain = combine(dfs)
    
    dfNorm, scaler = normalize(dfMain)

    markerList = list(set(dfNorm['marker']))
    
    print(markerList)

    # separate input, output, train and test

    '''
    for ident, df in dfs.items():
        
        df, scaler = normalize(df)
        input_indices, output_indices = split(df)

        in_df = np.array(df.iloc[input_indices[200][0]:input_indices[200][-1], :])
        out_df = np.array(df.iloc[output_indices[200][0]:output_indices[200][-1], :])
        back_df = scaler.inverse_transform(in_df)
    '''
# split across each copy separately
# manually train by epoch with a loop and set the fit epochs to 1
    # label by the label in above to reseparate when looping through; epoch then copy then batches
# batch_size should be the same as number of copies
# need to do by copy in order to use the stateful stuff


    ''' 
    df.to_csv("./entire.csv")
    in_df.to_csv("./inputExample.csv")
    out_df.to_csv("./outputExample.csv")
    
    values = df.values
    fig, axs = plt.subplots(len(df.columns))
    for i, ax in enumerate(axs):
        ax.plot(values[:,i])
        ax.set_title(list(df.columns)[i], loc='right')

    plt.savefig("data7.png")
    '''
