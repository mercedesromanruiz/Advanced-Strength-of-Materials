import sys
import math
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
    L = 2.0
    vertices = [
        (0.0, 0.0), # 0
        (L/2, 0.0), # 1
        (L, 0.0)] # 2

    cells = [
        (0, 1),
        (1, 2)]

    # vertexsets, with names and vertex labels
    # careful: if only one vertex, use a comma after the vertexlabel
    vertexsets = {}
    vertexsets["LPin"] = (0,)
    vertexsets["RPin"] = (2,)
    vertexsets["center"] = (1,)
    
    cellsets = {}
    cellsets["beams"] = (0, 1)


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
    E = 210e9
    A = 1e-3
    sigmae = 6e6

    elmtTypes["beams"] = FrameCreator({"E" : E,
                                       "A" : A,
                                       "I" : 1e-6,
                                       "W" : 1e-5,
                                       "sigmae": sigmae})

    # constrain a nodeset: (doflabel, value) for all nodes in the set
    constraints = defaultdict(list) # do not change this line
    constraints["LPin"].append((0, 0.0))
    constraints["LPin"].append((1, 0.0))
    constraints["RPin"].append((1, 0.0))
    # load a nodeset: (doflabel, value) for all nodes in the set
    loading = defaultdict(list) # do not change this line
    r = (A / math.pi)**(1/2)
    P = math.pi**3 / 2 * r**4 / L**2 *E
    loading["LPin"].append((0, P))
    loading["RPin"].append((0, -P))

    theModel = Model(theMesh, elmtTypes, constraints, loading)
    return theModel


def main():
    model = createModel()
    analysis = StaticLinearAnalysis(model)
    analysis.solve()
    model.printDetailed()



if __name__ == "__main__":
    main()
