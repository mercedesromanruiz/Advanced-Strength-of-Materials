#
# 2D truss element

import numpy as np
import math
from element import ElementCreator, Element


class TrussCreator(ElementCreator):
    """
    This is the class that contains the truss element template, that
    is, it contains the information that is shared by all the trusses
    that have the same properties.
    """

    def __init__(self, properties):
        """
        The constructor fo the TrussCreator is a dictionary of
        properties. The following properties are required:
            E: Young's modulus of the material
            A: cross section area.
            rho: density
            eigenstrain: from thermal, lack-of-fit, etc
        """
        self.YoungModulus = properties["E"]
        self.area = properties["A"]
        self.eigenstrain = properties.get("eigenstrain", 0.0)
        self.rho = properties.get("density", 1.0)
        self.g   = properties.get("gravity", 0.0)

        self.plottingShape = 1
        self.resultNames = ("epsilon", "sigma", \
                            "N", "stretch", "Vint", "buckling")


    def createElement(self, cell, nodes):
        return trussElement(cell, nodes, self)


    def getDOFsPerNode(self) -> int:
        return 2


    def print(self):
        print("Element type: truss")
        print("Young modulus:", self.YoungModulus)
        print("Cross section area: ", self.area)
        print("Eigenstrain:", self.eigenstrain)
        print("Density:", self.rho)
        print("Gravity:", self.g)


class trussElement(Element):

    def __init__(self, cell, nodes, eltype):
        """
        This is the class constructor. In this function, the information from the cell,
        the nodes, and the element type is manipulated and stored for future access. The
        variables that are created and can be accessed later on are:
           theNodes: the two nodes of the truss
           theType: the element type, with all the information provided by the user
           theCell: the cell (connecting the noces)
           L: the reference length of the bar
           volume: L * A

        All of these variables can be retrieved using self.theNodes, self.theType, etc.
        """
        self.theNodes = nodes
        self.theType = eltype
        self.theCell = cell

        nodei = self.theNodes[0]
        nodej = self.theNodes[1]
        r = nodej.coordinates() - nodei.coordinates()
        self.L = math.sqrt(r.dot(r))
        self.volume = self.L * eltype.area


    def Bmatrix(self):
        """
        Compute and return the element Bmatrix, the strain matrix.
        """
        nodei = self.theNodes[0]
        nodej = self.theNodes[1]

        r = nodej.coordinates() - nodei.coordinates()
        L = self.L
        d = r/L
        b = np.array([[-d[0]/L, -d[1]/L, d[0]/L, d[1]/L]])
        return b


    def bucklingCoeff(self, N):
        """
        arguments:
        N: axial force

        returns:
        Computes buckling coefficient phi for truss.
        phi < 1, bar is safe. The farther from 1, the better
        phi = 1, bar is buckling.
        phi > 1, bar has buckled!
        """
        phi = 0.0
        return phi


    def integrateEnergy(self):
        """
        Integrates the total potential energy in the truss element.
        Returns it value.
        """
        E = self.theType.YoungModulus
        A = self.theType.area
        e0 = self.theType.eigenstrain
        rho = self.theType.rho
        g = self.theType.g

        nodei = self.theNodes[0]
        nodej = self.theNodes[1]
        b = self.Bmatrix()

        disp = np.array([nodei.U[0], nodei.U[1], nodej.U[0], nodej.U[1]])
        epsilon = np.matmul(b, disp)[0]
        weight = rho * g * A * self.L
        fext = np.array([0.0, -0.5*weight, 0.0, -0.5*weight])

        energy = 0.5 * E * A * (epsilon - e0) * (epsilon - e0) * self.L - np.dot(disp,fext)
        return energy


    def integrateDEnergy(self):
        """
        Compute and return the contribution of the element to
        the equilibrium equations.
        """
        E = self.theType.YoungModulus
        A = self.theType.area
        e0 = self.theType.eigenstrain
        rho = self.theType.rho
        g = self.theType.g

        nodei = self.theNodes[0]
        nodej = self.theNodes[1]
        b = self.Bmatrix()

        disp = np.array([nodei.U[0], nodei.U[1], nodej.U[0], nodej.U[1]])
        epsilon = np.matmul(b, disp)[0]
        N = E * A * (epsilon - e0)

        BTs = b.T * N * self.L
        weight = rho * g * A * self.L
        fext = np.array([[0.0, -0.5*weight, 0.0, -0.5*weight]])

        return BTs - fext.T


    def integrateDDEnergy(self):
        """
        Compute and return the element stiffness matrix.
        """
        E = self.theType.YoungModulus
        A = self.theType.area
        e0 = self.theType.eigenstrain

        b = self.Bmatrix()
        K = E * A * np.outer(b,b) * self.L

        return K


    def print(self):
        """
        Print the information about the element.
        """
        self.theType.print()

        nodei = self.theNodes[0]
        nodej = self.theNodes[1]
        b = self.Bmatrix()

        disp = np.array([[nodei.U[0], nodei.U[1], nodej.U[0], nodej.U[1]]]).T
        epsilon = np.matmul(b, disp)[0,0]

        E = self.theType.YoungModulus
        A = self.theType.area
        L = self.L
        N = E * A * epsilon

        Vint = 0.5 * E * A * epsilon * epsilon * self.L
        buck = self.bucklingCoeff(N)

        print("Node 1 DOFS:", nodei.DOFS)
        print("Node 2 DOFS:", nodej.DOFS)

        print("epsilon: ", epsilon, ", sigma: ", E*epsilon)
        print("N: ", E * A * epsilon, "stretch: ", epsilon*L)
        print("internal energy: ", Vint)
        print("Buckling criterion: ", buck)


    def result(self, name):
        """
        Computes a scalar called 'name' that depends
        on the state of the truss. It is used by the plotting functions.
        """
        nodei = self.theNodes[0]
        nodej = self.theNodes[1]
        L = self.L
        b = self.Bmatrix()
        disp = np.array([[nodei.U[0], nodei.U[1], nodej.U[0], nodej.U[1]]]).T
        epsilon = np.matmul(b, disp)[0,0]
        E = self.theType.YoungModulus
        A = self.theType.area
        N = E * A * epsilon

        buck = self.bucklingCoeff(N)

        if (name == "epsilon"):
            r = epsilon

        elif (name == "sigma"):
            r = E * epsilon

        elif (name == "N"):
            r = E * A * epsilon

        elif (name == "stretch"):
            r = epsilon * L

        elif (name == "Vint"):
            r = 0.5 * E * A * epsilon * epsilon * self.L

        elif (name == "buckling"):
            r = buck

        else:
            r = 0.0

        return r
