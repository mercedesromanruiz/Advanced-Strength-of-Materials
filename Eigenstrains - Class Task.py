# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 13:21:23 2024

@author: merce
"""

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
    
    vertices = [
        (0.0, 0.0),
        (15.0, 0.0)]

    cells = [
        (0, 1)]

    # vsets, with names and vertex labels
    # careful: if only one vertex, use a comma after the vertexlabel
    vsets = {}
    vsets["Izq"] = (0,)
    vsets["Der"] = (1,)

    # csets, with names and cell labels
    # careful: if only one cell, use a comma after the cell label
    cellsets = {}
    cellsets["bars"] = (0,)


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
    
    elmtTypes = {} # do not change this line
    elmtTypes["bars"] = TrussCreator({"E":210.0e9,
                                      "A":0.02,
                                      "eigenstrain":0.03})

    # constrain a vertexset: (doflabel, value) for all nodes in the set
    constraints = defaultdict(list) # do not change this line
    constraints["Izq"].append((0, 0.0))
    constraints["Izq"].append((1, 0.0))
    
    constraints["Der"].append((0, 0.0))
    constraints["Der"].append((1, 0.0))

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
