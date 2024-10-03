import sys
sys.path.append('src')

import math
from fem import Mesh, Model, StaticLinearAnalysis
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
    #           mesh geometry: vertices, cells and nodesets
    #--------------------------------------------------------
    a = 2
    vertices = [
        (0.0, 0.0),
        (  a, 0.0),
        (a/2, (a/2)*math.tan(math.radians(30)))]

    cells = [
        (0, 2),
        (1, 2)]

    # vsets, with names and vertex labels
    # careful: if only one vertex, use a comma after the vertexlabel
    vsets = {}
    vsets["support"] = (0, 1)
    vsets["tip"] = (2,)

    # csets, with names and cell labels
    # careful: if only one cell, use a comma after the cell label
    cellsets = {}
    cellsets["mybars"] = (0, 1)


    # --------------------------------------------------------
    #               Create global data structures
    # --------------------------------------------------------
    theMesh = Mesh(vertices, cells, vsets, cellsets)
    theMesh.print()

    # --------------------------------------------------------
    #  Model: elements types, constraints, and loading on the mesh
    # --------------------------------------------------------

    # define and assign the element types
    # 'E' is Young's modulus of the material
    # 'A' is the cross section area
    EA = 2e6
    elmtTypes = {} # do not change this line
    elmtTypes["mybars"] = TrussCreator({"E" : EA, "A": 1.0})

    # constrain a vertexset: (doflabel, value) for all nodes in the set
    constraints = defaultdict(list) # do not change this line
    constraints["support"].append((0, 0.0))
    constraints["support"].append((1, 0.0))

    # load a vertexset: (doflabel, value) for all nodes in the set
    F = 1000
    loading = defaultdict(list) # do not change this line
    loading["tip"].append((0, F))

    theModel = Model(theMesh, elmtTypes, constraints, loading)
    return theModel


def main():
    model = createModel()
    model.print()
    analysis = StaticLinearAnalysis(model)
    analysis.solve()
    model.printDetailed()


if __name__ == "__main__":
    main()
