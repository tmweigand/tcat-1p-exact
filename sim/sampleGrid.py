import numpy as np

L = [[1.1,5],[-2.5,2.5]]

nX = 20
nY = 20

dX = (L[0][1]-L[0][0])/nX
dY = (L[1][1]-L[1][0])/nY

for i in range(0,nX+1):
    for j in range(0,nY+1):
        print("(",L[0][0]+i*dX,L[1][0]+j*dY,0.,")")