import os
import numpy as np


def sortFunc(d):
    return int(d.split("/")[2])

rootdir = 'postProcessing/data'
dataFile = []
for file in os.listdir(rootdir):
    d = os.path.join(rootdir, file)
    if os.path.isdir(d):
        dataFile.append(d)
    
dataFile.sort(key=sortFunc)
numData = len(dataFile)

data = {}
for n,dFile in enumerate(dataFile[0:10]):
    arr = np.loadtxt(dFile+"/data_U.csv",delimiter=",",skiprows=1)
    ID = sortFunc(dFile)
    data[n] = {'data':arr,'time':ID}
