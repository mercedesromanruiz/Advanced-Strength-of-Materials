import sys
sys.path.append('src')

import math
from fem import Mesh, Model, StaticLinearAnalysis
from collections import defaultdict
from elements.truss import TrussCreator
from elements.frame import FrameCreator


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
    vertices = [
        (0.0, 15.0),
        (30.0, 15.0),
        (15.0, 20.0),
        (12.5, 15.0),
        (17.5, 15.0),
        (7.5, 17.5),
        (22.5, 17.5)]

    cells = [
        (0, 3),
        (0, 5),
        (3, 4),
        (3, 2),
        (3, 5),
        (4, 1),
        (4, 2),
        (4, 6),
        (2, 6),
        (5, 2),
        (6, 1)]

    # vsets, with names and vertex labels
    # careful: if only one vertex, use a comma after the vertexlabel
    vertexsets = {}
    vertexsets["support_fixed"] = (0,)
    vertexsets["support_mobile"] = (1,)
    vertexsets["center_left"] = (3,)
    vertexsets["center_right"] = (4,)

    # csets, with names and cell labels
    # careful: if only one cell, use a comma after the cell label
    cellsets = {}
    cellsets["deck"] = (0, 2, 5)
    cellsets["diagonals"] = (1, 3, 4, 6, 7, 8, 9, 10)

    # --------------------------------------------------------
    #               Create global data structures
    # --------------------------------------------------------
    theMesh = Mesh(vertices, cells, vertexsets, cellsets)
    #theMesh.print()

    # --------------------------------------------------------
    #  Model: elements types, constraints, and loading on the mesh
    # --------------------------------------------------------

    # define and assign the element types
    # 'E' is Young's modulus of the material
    # 'A' is the cross section area

    elmtTypes = {} # do not change this line

    # Steel
    sigmae = 180e6 # Strength
    density = 7800 # Density
    EA = 2.1e11 # Stiffness

    A = 1e-2
    r = (A / math.pi)**(1/2)
    I = math.pi / 4 * r**4
    W = I / r
    
    
    elmtTypes["deck"] = FrameCreator({"E" : EA / A,
                                      "A" : A,
                                      "I" : I,
                                      "fx" : 10,
                                      "fy": -100,
                                      "gravity" : 9.8,
                                      "density": density,
                                      "sigmae": sigmae,
                                      "W" : W})
    
    elmtTypes["diagonals"] = FrameCreator({"E" : EA / A,
                                           "A" : A,
                                           "I" : I,
                                           "fx": 10,
                                           "fy": 0,
                                           "gravity" : 9.8,
                                           "density": density,
                                           "sigmae": sigmae,
                                           "W" : W})


    # constrain a vertexset: (doflabel, value) for all nodes in the set
    constraints = defaultdict(list) # do not change this line
    constraints["support_fixed"].append((0, 0.0))
    constraints["support_fixed"].append((1, 0.0))
    constraints["support_mobile"].append((1, 0.0))

    # load a vertexset: (doflabel, value) for all nodes in the set
    loading = defaultdict(list) # do not change this line

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
