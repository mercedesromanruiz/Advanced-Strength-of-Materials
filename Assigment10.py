import sys
sys.path.append('src')


import numpy as np
from fem import Mesh, Model, StaticLinearAnalysis
from collections import defaultdict
from elements.frame import FrameCreator
from elements.pointmass import PointmassCreator

def printModelInfo(model):
    print("Model volume: ", model.integrateVolume())
    maxdofs = model.maximumNodalDOFs()
    print("Maximum vertical displacement:", maxdofs[1])


def elementWithMinEnergy(model):
    minEner = 1000000
    thin = model.eltypes["thin beams"]

    for e, el in model.elements.items():
        elener = el.integrateEnergy(1.0)
        if elener < minEner and el.theType != thin:
            minEner = elener
            worstElem = el

    return worstElem


# create a grid of vertices
lx = 20
ly = 10
nx = 20
ny = 15
dx = lx/(nx-1)
dy = ly/(ny-1)

vertices = []
for j in range(ny):
    for i in range(nx):
        vertices.append((i*dx, j*dy))

# create a mesh of beams
cells = []
for i in range(nx-1):
    for j in range(ny-1):
        cells.append((j*nx+i, j*nx+i+1))
        cells.append((j*nx+i, (j+1)*nx+i))
        cells.append((j*nx+i, (j+1)*nx+i+1))
        cells.append((j*nx+i+1, (j+1)*nx+i))
    cells.append(((ny-1)*nx+i, (ny-1)*nx+i+1))

for j in range(ny-1):
    cells.append(((j+1)*nx-1, (j+2)*nx-1))

# vertexsets, with names and vertex labels
# careful: if only one vertex, use a comma after the vertexlabel
# select suppports and nodes for loading
ls = []
rs = []
ll = []
for nv, v in enumerate(vertices):
    if (v[0] == 0):
       ls.append(nv)
    if (v[0] == 20):
        rs.append(nv)
    if (v[1] == ly and v[0] >= 7 and v[0]<=13):
        ll.append(nv)

vertexsets = {}
vertexsets["left support"] = tuple(ls)
vertexsets["right support"] = tuple(rs)
vertexsets["for loading"] = tuple(ll)

cellsets = {}
cellsets["thick beams"] = range(len(cells))
cellsets["thin beams"] = ()


# --------------------------------------------------------
#               Create global data structures
# --------------------------------------------------------
theMesh = Mesh(vertices, cells, vertexsets, cellsets)
theMesh.print()


# --------------------------------------------------------
#  Model: elements types, constraints, and loading on the mesh
# --------------------------------------------------------

# assign one element type to each cellset
elmtTypes = {}
elmtTypes["thick beams"] = FrameCreator({"E" : 30e9,
                                   "A" : 7.069e-4,
                                   "I" : 3.976e-8,
                                   "density" : 2000,
                                   "gravity" : 9.8,
                                   "W" : 2.651e-6,
                                   "yield": 60e6})

elmtTypes["thin beams"] = FrameCreator({"E" : 30e9,
                                   "A" : 7.069e-8,
                                   "I" : 3.976e-16,
                                   "density" : 2000,
                                   "gravity" : 9.8,
                                   "W" : 2.651e-6,
                                   "yield": 60e6})


theModel = Model(theMesh, elmtTypes)


theModel.addConstraint(vertexset="left support" , dof=0, value=0.0)
theModel.addConstraint(vertexset="left support" , dof=1, value=0.0)
theModel.addConstraint(vertexset="left support" , dof=2, value=0.0)

theModel.addLoading(vertexset="right support", dof=1, value=-60000.0)


analysis = StaticLinearAnalysis(theModel)
theModel.print()
analysis.solve()
printModelInfo(theModel)
vol0 = theModel.integrateVolume()
thin = theModel.eltypes["thin beams"]
maxvertical = theModel.maximumNodalDOFs()[1]
print("\n*** Max vertical displacement", maxvertical)

saving = 0
count = 2
while saving < 25:
    v = theModel.integrateVolume()
    saving = (vol0-v)/vol0*100
    print("*** Volume saving (%):", saving)
    el = elementWithMinEnergy(theModel)
    el.theType = thin
    analysis.resetSolution()
    analysis.solve(label=count)
    maxvertical = theModel.maximumNodalDOFs()[1]
    print("\n*** Max vertical displacement", maxvertical)
    count += 1

print("\n\nOptimal structure found")
print(">>>> Max vertical displacement", maxvertical)
print(">>>> Percentage volume saving: ", saving)
print(">>>> Volume", theModel.integrateVolume())
