import os
import numpy as np

folder = "postProcessing"
import os
 
rootdir = 'postProcessing'
for file in os.listdir(rootdir):
    d = os.path.join(rootdir, file)
    if os.path.isdir(d):
        print(d)