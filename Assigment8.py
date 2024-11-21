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
    (0.0, 0.0), # 0 # Left Base
    (0.0, L/3),
    (0.0, 2*L/3),
    (0.0, L),
    (L/3, L),
    (2*L/3, L),
    (L, L), # 6 # Top Right
    (L, 2*L/3),
    (L, L/3),
    (L, 0)] # 9 # Right Base

cells = [
    (0, 1),
    (1, 2),
    (2, 3),
    (3, 4),
    (4, 5),
    (5, 6),
    (6, 7),
    (7, 8),
    (8, 9),
    ]

# vertexsets, with names and vertex labels
# careful: if only one vertex, use a comma after the vertexlabel
vertexsets = {}
vertexsets["Clamp"] = (0, 9)
vertexsets["TopRight"] = (6,)

cellsets = {}
cellsets["SteelFrame"] = (0, 1, 2, 3, 4, 5, 6, 7, 8)


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
side = 8e-2
sigmae = 180e6 # Strength
density = 7800 # Density
EA = 2.1e11 # Stiffness
A = side * side
I = 1 / 12 * side**4

elmtTypes["SteelFrame"] = FrameCreator({"E" : EA / A,
                                        "A" : A,
                                        "I" : I,
                                        "gravity" : 9.8,
                                        "density": density,
                                        "sigmae": sigmae})

# constrain a nodeset: (doflabel, value) for all nodes in the set
constraints = defaultdict(list) # do not change this line
constraints["Clamp"].append((0, 0.0))
constraints["Clamp"].append((1, 0.0))
constraints["Clamp"].append((2, 0.0))


loading = {}
F = 1000
loading["TopRight"] = Load(dof=1, value=F, scaling='t*(t<=10)')

theModel = Model(theMesh, elmtTypes, constraints, loading)

dt = 0.01
tf = 10.0
analysis = TransientAnalysis(theModel, dt, tf)
analysis.solve()
#model.printDetailed()