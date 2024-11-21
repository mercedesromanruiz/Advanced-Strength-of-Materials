import sys
sys.path.append('src')
import math

from fem import *
from collections import defaultdict
from elements.frame import FrameCreator

#--------------------------------------------------------
#           mesh geometry: vertices, cells and nodesets
#--------------------------------------------------------
L = 2.0
vertices = [
    (0.0, 0.0),
    (L/2, 0.0),
    (L, 0.0)]

cells = [(0, 1), (1,2)]

# vertexsets, with names and vertex labels
# careful: if only one vertex, use a comma after the vertexlabel
vertexsets = {}
vertexsets["clamp"] = (0,)
vertexsets["tip"] = (2,)

cellsets = {}
cellsets["beams"] = (0,1)


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
elmtTypes["beams"] = FrameCreator({"E" : 210e9,
                                   "A" : 1e-2,
                                   "I" : 1e-6,
                                   "W" : 1e-5,
                                   "sigmae": 2e6,
                                   "density" : 400000.0,
                                   "gravity" : 0.0})

# constrain a nodeset: (doflabel, value) for all nodes in the set
constraints = defaultdict(list) # do not change this line
constraints["clamp"].append((0, 0.0))
constraints["clamp"].append((1, 0.0))
constraints["clamp"].append((2, 0.0))


loading = {}
loading["tip"] = Load(dof=1, value=-20000.0, scaling='t*(t<=2)')

theModel = Model(theMesh, elmtTypes, constraints, loading)

dt = 0.01
tf = 10.0
analysis = TransientAnalysis(theModel, dt, tf)
analysis.solve()
#model.printDetailed()
