import meshio
import numpy as np
import matplotlib.pyplot as plt
from math import *

def procMsh(file_name):
    
    msh = meshio.read(file_name)

    X = msh.points[:,0]
    Y = msh.points[:,1]

    IEN = msh.cells[1].data
    IENbound = msh.cells[0].data

    IENboundTypeElem = list(msh.cell_data['gmsh:physical'][0])
    
    #print(msh.cell_data['gmsh:physical'])
    
    #print(IENboundTypeElem)

    boundNames= list(msh.field_data.keys())

    #IENboundElem = [boundNames[elem] for elem in IENboundTypeElem]

    IENboundElemNames = []

    for i in IENboundTypeElem:
        if i == 135:
            IENboundElemNames.append(boundNames[0])
        if i == 136:
            IENboundElemNames.append(boundNames[1])
        if i == 137:
            IENboundElemNames.append(boundNames[2])
        if i == 138:
            IENboundElemNames.append(boundNames[3])
        if i == 139:
            IENboundElemNames.append(boundNames[4])
        if i == 140:
            IENboundElemNames.append(boundNames[5])

    return [X, Y, IEN, IENbound, IENboundElemNames]

if __name__ == '__main__':
    
    procMsh('mshTri2.msh')
    
    #print(procMsh('malha.msh')[3])
