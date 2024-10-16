import sys
sys.path.append('../../../src')

from fem import Mesh, Model
from collections import defaultdict
from elements.truss import TrussCreator


def createModel():
    """
    The goal of one example is to create a Model. This model will have
    to include the Mesh (with its components - see below), the type of
    elements (the equations that represent its behaviour), and the
    actions on the model (boundary conditions and loading).
    """

    #--------------------------------------------------------
    #           mesh geometry: vertices, cells and sets
    #--------------------------------------------------------
    vertices = [
        (0.0, 0.0),
        (15.0, 0.0),
        (30.0, 0.0),
        (15.0, 15.0)]

    cells = [
        (0, 1),
        (0, 3),
        (1, 3),
        (1, 2),
        (3, 2)]

    # vertexsets, with names and vertex labels
    # careful: if only one vertex, use a comma after the vertexlabel
    vertexsets = {}
    vertexsets["support"] = (0, 2)
    vertexsets["center"] = (1,)

    cellsets = {}
    cellsets["bars"] = (0, 1, 2, 3, 4)

    # --------------------------------------------------------
    #               Create global data structures
    # --------------------------------------------------------
    theMesh = Mesh(vertices, cells, vertexsets, cellsets)
    #theMesh.print()


    # --------------------------------------------------------
    #  Model: elements types, constraints, and loading on the mesh
    # --------------------------------------------------------

    # assign one element type to each cellset
    elmtTypes = {}
    elmtTypes["bars"] = TrussCreator({"E": 210.0e9,
                                      "A": 1e-4,
                                      "density": 1700,
                                      "gravity": 9.8})

    # constrain a vertexset: (doflabel, value) for all nodes in the set
    constraints = defaultdict(list) # do not change this line
    constraints["support"].append((0, 0.0))
    constraints["support"].append((1, 0.0))

    # load a vertexset: (doflabel, value) for all nodes in the set
    loading = defaultdict(list) # do not change this line
    loading["center"].append((1, -10000.0))

    theModel = Model(theMesh, elmtTypes, constraints, loading)
    return theModel
