import sys
sys.path.append('src')

from fem import Mesh, Model, StaticLinearAnalysis
from collections import defaultdict
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
        (0.0, 0.0), # 0 # Left Base
        (0.0, 2.44/2), # 1
        (0.0, 2.44), # 2
        (7.32/2, 2.44), # 3 # Center Top Beam -> Loading
        (7.32,   2.44), # 4
        (7.32, 2.44/2), # 5
        (7.32, 0.0)] # 6 # Right Base

    cells = [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 6)]

    # vertexsets, with names and vertex labels
    # careful: if only one vertex, use a comma after the vertexlabel
    vertexsets = {}
    vertexsets["Clamp"] = (0, 6)
    vertexsets["CenterTop"] = (3,)
    
    cellsets = {}
    cellsets["WoodenBeams"] = (0, 1, 2, 3, 4, 5)


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
    elmtTypes["WoodenBeams"] = FrameCreator({"E" : 12e9,
                                       "A" : 0.12 * 0.12,
                                       "I" : 1/12 * 0.12 * 0.12**3})

    # constrain a nodeset: (doflabel, value) for all nodes in the set
    constraints = defaultdict(list) # do not change this line
    constraints["Clamp"].append((0, 0.0))
    constraints["Clamp"].append((1, 0.0))
    constraints["Clamp"].append((2, 0.0))

    # load a nodeset: (doflabel, value) for all nodes in the set
    loading = defaultdict(list) # do not change this line
    loading["CenterTop"].append((1, -80 * 9.8))

    theModel = Model(theMesh, elmtTypes, constraints, loading)
    return theModel

def main():
    model = createModel()
    analysis = StaticLinearAnalysis(model)
    analysis.solve()
    model.printDetailed()


if __name__ == "__main__":
    main()
