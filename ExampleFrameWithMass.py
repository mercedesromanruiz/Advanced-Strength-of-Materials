import sys
sys.path.append('src')

from fem import Mesh, Model, StaticLinearAnalysis
from collections import defaultdict
from elements.frame import FrameCreator
from elements.truss import TrussCreator
from elements.pointmass import PointmassCreator

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
    L = 2.0
    vertices = [
        (0.0, 0.0), # 0
        (0.5, 0.0), # 1
        (1.0, 0.0), # 2
        (0.0, 0.8), # 3
        (0.5, 0.8), # 4
        (1.0, 0.8), # 5
        (0.0, 1.3), # 6
        (0.5, 1.3), # 7
        (1.0, 1.3)] # 8

    cells = [
        (0, 3), # 0
        (1, 4), # 1
        (2, 5), # 2
        (3, 6), # 3
        (4, 7), # 4
        (5, 8), # 5
        (3, 4), # 6
        (4, 5), # 7
        (6, 7), # 8
        (7, 8), # 9
        (3, 7), # 10
        (5, 7), # 11
        (3,),   # 12
        (4,),   # 13
        (5,),   # 14
        (8,)    # 15
        ]

    # vertexsets, with names and vertex labels
    # careful: if only one vertex, use a comma after the vertexlabel
    vertexsets = {}
    vertexsets["floor"] = (0, 1, 2)
    vertexsets["top-right-corner"] = (8,)

    cellsets = {}
    cellsets["columns"] = (0, 1, 2, 3, 4, 5)
    cellsets["horizontal"] = (6, 7, 8, 9)
    cellsets["diagonals"] = (10, 11)
    cellsets["smallmasses"] = (12, 13, 14)
    cellsets["bigmass"] = (15,)


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
    elmtTypes["columns"] = FrameCreator({"E" : 210e9,
                                       "A" : 1e-2,
                                       "I" : 3.3e-6,
                                       "density" : 1.0,
                                       "gravity" : 9.8})
    elmtTypes["horizontal"] = FrameCreator({"E" : 210e9,
                                       "A" : 5e-3,
                                       "I" : 3.3e-5,
                                       "density" : 1.0,
                                       "gravity" : 9.8})
    elmtTypes["diagonals"] = TrussCreator({"E" : 210e9,
                                       "A" : 5e-3,
                                       "I" : 3.3e-5,
                                       "density" : 1.0,
                                       "gravity" : 9.8})
    elmtTypes["smallmasses"] = PointmassCreator({"mass" : 300.0,
                                           "gravity" : 9.8})
    elmtTypes["bigmass"] = PointmassCreator({"mass" : 600.0,
                                           "gravity" : 9.8})

    theModel = Model(theMesh, elmtTypes)
    theModel.addConstraint(vertexset = "floor", dof=0, value=0.0)
    theModel.addConstraint(vertexset = "floor", dof=1, value=0.0)
    theModel.addConstraint(vertexset = "floor", dof=2, value=0.0)
    theModel.addLoading(vertexset="top-right-corner", dof=0, value=-10000.0)
    return theModel


def main():
    model = createModel()
    analysis = StaticLinearAnalysis(model)
    analysis.solve()
    model.printDetailed()


if __name__ == "__main__":
    main()
