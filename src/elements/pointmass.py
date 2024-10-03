#
# 2D pointmass element

import numpy as np
import math
from element import ElementCreator, Element


class PointmassCreator(ElementCreator):
    """This is the class that contains the pointmass element template, that
    is, it contains the information that is shared by all the pointmasses
    that have the same properties.
    """

    def __init__(self, properties):
        """
        The constructor fo the TrussCreator is a dictionary of
        properties. The following properties are required:
            mass: concentrated mass
            gravity: acceleration
        """
        self.m = properties.get("mass", 0.0)
        self.g = properties.get("gravity", 0.0)

        self.plottingShape = 0
        self.resultNames = ()


    def createElement(self, cell, nodes):
        return pointmassElement(cell, nodes, self)


    def getDOFsPerNode(self) -> int:
        return 2


    def print(self):
        print("Element type: pointmass")
        print("Concentrated mass:", self.m)
        print("Gravity: ", self.g)


class pointmassElement(Element):

    def __init__(self, cell, nodes, eltype):
        """
        This is the class constructor. In this function, the information from the cell,
        the nodes, and the element type is manipulated and stored for future access. The
        variables that are created and can be accessed later on are:
           theNodes: the two nodes of the pointmass
           theType: the element type, with all the information provided by the user
           theCell: the cell (connecting the noces)
           L: the reference length of the bar
           volume: L * A

        All of these variables can be retrieved using self.theNodes, self.theType, etc.
        """
        self.theNodes = nodes
        self.theType = eltype
        self.theCell = cell
        self.volume = 0.0


    def integrateEnergy(self):
        """
        Integrates the total potential energy in the pointmass element.
        Returns it value.
        """
        m = self.theType.m
        g = self.theType.g

        node = self.theNodes[0]
        disp = np.array([node.U[0], node.U[1]])
        weight = m * g
        fext = np.array([0.0, -weight])

        energy = -np.dot(disp,fext)
        return energy


    def integrateDEnergy(self):
        """
        Compute and return the contribution of the element to
        the equilibrium equations.
        """
        m = self.theType.m
        g = self.theType.g

        node = self.theNodes[0]
        disp = np.array([node.U[0], node.U[1]])
        weight = m * g
        fext = np.array([0.0, -weight])
        return -fext.T


    def integrateDDEnergy(self):
        """
        Compute and return the element stiffness matrix.
        """
        K = np.zeros((2,2))
        return K


    def print(self):
        """
        Print the information about the element.
        """
        self.theType.print()

        node = self.theNodes[0]
        disp = np.array([[node.U[0], node.U[1]]]).T
        print("Displacement: Ux:", disp[0], ", Uy:", disp[1])


    def result(self, name):
        """
        Return results to be plotted with paraview. None for
        the pointmass
        """
        return 0.0
