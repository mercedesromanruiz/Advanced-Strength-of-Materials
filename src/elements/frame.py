#
# 2D frame element

import numpy as np
import math
from element import ElementCreator, Element
from element import quadraturePoints, evaluateShapeFunctions, GaussLobattoPoints


class FrameCreator(ElementCreator):
    """This is the class that contains the frame element template, that
    is, it contains the information that is shared by all the frames
    that have the same properties.
    """

    def __init__(self, properties):
        """
        The constructor fo the FrameCreator is a dictionary of
        properties. The following properties are required:
            E: Young's modulus of the material
            A: cross section area.
            I: inertia
            W: section modulus
            sigmae: elastic limit
            density: density (per unit volume)
            gravity: vertical and downwards
            eigenstrain: from thermal, lack-of-fit, etc
            safetybuck: admissible safety criterion for buckling
            safetysigma: admissible safety criterion for stress

        """
        self.YoungModulus = properties["E"]
        self.area = properties["A"]
        self.inertia = properties["I"]
        self.modulus = properties.get("W", 0.0)
        self.sigmae  = properties.get("sigmae", -1.0)
        self.rho     = properties.get("density", 1.0)
        self.g       = properties.get("gravity", 0.0)
        self.eigenstrain = properties.get("eigenstrain", 0.0)
        self.safetybuck = properties.get("safetybuck", -1.0)
        self.safetysigma = properties.get("safetysigma", -1.0)

        # these are the variables that will be plotted by paraview.
        # all of them will have a single value for each beam, so one
        # should not expect smooth colors
        self.resultNames = ("N", "epsilon",
                            "M", "kappa",
                            "safety")


    def createElement(self, cell, nodes):
        return frameElement(cell, nodes, self)


    def getDOFsPerNode(self) -> int:
        """
        Each node of the plane frame has three dofs: (ux, uy, theta)
        """
        return 3

    def print(self):
        print("Element type: frame")
        print("Young modulus:", self.YoungModulus)
        print("Cross section area:", self.area)
        print("Cross section inertia:", self.inertia)
        print("Section modulus:", self.modulus)
        print("Density:", self.rho)
        print("Gravity:", self.g)

        if self.eigenstrain > 0:
            print("Eigenstrain:", self.eigenstrain)

        printIfPositive("Elastic limit: ", self.sigmae)
        printIfPositive("Admissible safety criterion for buckling: ", self.safetybuck)
        printIfPositive("Admissible safety criterion for stresses: ", self.safetysigma)


class frameElement(Element):

    def __init__(self, cell, nodes, eltype):
        """
        This is the class constructor. In this function, the information from the cell,
        the nodes, and the element type is manipulated and stored for future access. The
        variables that are created and can be accessed later on are:
           theNodes: the two nodes of the frame
           theType: the element type, with all the information provided by the user
           theCell: the cell (connecting the noces)
           L: the reference length of the bar
           volume: L * A

        All of these variables can be retrieved using self.theNodes, self.theType, etc.
        """
        self.theNodes = nodes
        self.theType = eltype
        self.theCell = cell

        node0 = self.theNodes[0]
        node1 = self.theNodes[1]
        r = node1.coordinates() - node0.coordinates()
        self.L = math.sqrt(r.dot(r))
        self.volume = self.L * eltype.area
        self.alpha = math.atan2(r[1], r[0])


    def Bmatrix(self, x):
        """
        Compute and return the element Bmatrix, the strain matrix at x.
        """
        node0 = self.theNodes[0]
        node1 = self.theNodes[1]
        r = node1.coordinates() - node0.coordinates()
        L = self.L
        xi = x/L

        bh = np.array([[-1.0/L, 0.0, 0.0, 1.0/L, 0.0, 0.0],
                       [0.0, L*L*(12*xi-6.0), 2*L**3*(3*xi-2),
                        0.0, L**2*(6-12*xi), 2*L**3*(3*xi-1)]])

        r = self.rotation()
        b = np.matmul(bh, r)
        return b


    def integrateEnergy(self):
        """
        Integrates the total potential energy in the frame element.
        Returns it value.
        """
        k = self.localStiffness()
        r = self.rotation()

        node0 = self.theNodes[0]
        node1 = self.theNodes[1]

        U = np.array([node0.U[0], node0.U[1], node0.U[2],
                      node1.U[0], node1.U[1], node1.U[2]])

        bigK = np.matmul(r.T, np.matmul(k,r))

        f = self.localForce()
        bigF = np.matmul(r.T, f)

        energy = 0.5*np.matmul(U, np.matmul(bigK, U)) - np.dot(bigF,U)
        return energy


    def integrateDEnergy(self):
        """
        Compute and return the contribution of the element to
        the equilibrium equations.
        """
        k = self.localStiffness()
        r = self.rotation()

        node0 = self.theNodes[0]
        node1 = self.theNodes[1]

        U = np.array([node0.U[0], node0.U[1], node0.U[2],
                      node1.U[0], node1.U[1], node1.U[2]])

        bigK = np.matmul(r.T, np.matmul(k,r))

        f = self.localForce()
        bigF = np.matmul(r.T, f)
        res = np.matmul(bigK, U) - bigF
        return res


    def integrateDDEnergy(self):
        """
        Compute and return the element stiffness matrix.
        """
        k = self.localStiffness()
        r = self.rotation()
        bigK = np.matmul(r.T, np.matmul(k,r))

        return bigK


    def integrateVolume(self):
        node0 = self.theNodes[0]
        node1 = self.theNodes[1]
        r = node1.coordinates() - node0.coordinates()
        L = math.sqrt(r.dot(r))
        vol = L * self.theType.area
        return vol


    def localForce(self):
        rho = self.theType.rho
        g = self.theType.g
        L = self.L
        f = rho*g*math.sin(self.alpha)
        q = rho*g*math.cos(self.alpha)

        fv = np.array([-f*L/2, -q*L/2, -q*L*L/12, -f*L/2, -q*L/2, q*L*L/12])
        return fv.T


    def localStiffness(self):
        E = self.theType.YoungModulus
        A = self.theType.area
        II = self.theType.inertia
        L = self.L
        L2 = L*L
        L3 = L2*L
        EA = E*A
        EI = E*II

        k = np.array(
           [[ EA/L,         0,        0, -EA/L,         0,       0],
            [    0,  12*EI/L3,  6*EI/L2,     0, -12*EI/L3, 6*EI/L2],
            [    0,   6*EI/L2,   4*EI/L,     0,  -6*EI/L2,  2*EI/L],
            [-EA/L,         0,        0,  EA/L,         0,       0],
            [    0, -12*EI/L3, -6*EI/L2,     0,  12*EI/L3, -6*EI/L2],
            [    0,   6*EI/L2,   2*EI/L,     0,  -6*EI/L2, 4*EI/L]])
        return k


    def print(self):
        """
        Print the state at the mid-span of the element.
        """
        self.theType.print()
        E = self.theType.YoungModulus
        A = self.theType.area
        II = self.theType.inertia
        L = self.L
        W = self.theType.modulus

        node0 = self.theNodes[0]
        node1 = self.theNodes[1]

        U = np.array([node0.U[0], node0.U[1], node0.U[2],
                      node1.U[0], node1.U[1], node1.U[2]])

        # calculate quantities at the left, center and right
        quadPoints = GaussLobattoPoints(3)
        sigmamax = float('-inf')
        sigmamin = float('inf')
        Gstress = float('inf')
        Gbuck = float('inf')

        for k, qp in enumerate(quadPoints):
            x = 0.5*(qp.x + 1)*L
            B = self.Bmatrix(x)

            [epsilon, kappa] = np.matmul(B, U)
            N = epsilon * E * A
            M = kappa * E * II

            Gstress_k, Gbuck_k = self.safetyFactors(N, M)
            Gstress = min(Gstress, Gstress_k)
            Gbuck = min(Gbuck, Gbuck_k)

            # maximum stress. This needs to be programmed
            sigmamin_k, sigmamax_k = self.stressMinMax(N, M)
            sigmamax = max(sigmamax, sigmamax_k)
            sigmamin = min(sigmamin, sigmamin_k)

            print("\n\tSection", k, " in beam at x:", x)
            print("\tAxial force    : {:6g}".format(N))
            print("\tAxial strain   : {:6g}".format(epsilon))
            print("\tBending moment : {:6g}".format(M))
            print("\tCurvature      : {:6g}".format(kappa))
            print("\tMaximum stress : {:6g}".format(sigmamax_k))
            print("\tMinimum stress : {:6g}".format(sigmamin_k))
            print("\tStress safety  : {:6g}".format(Gstress_k))
            print("\tBuckling safety: {:6g}".format(Gbuck_k))

        # buckling safefy factor for the whole beam
        print("\nBuckling safety factor: {0:.6g}".format(Gbuck))

        # stress safefy factor for the whole beam
        print("Stress safety factor: {0:.6g}".format(Gstress))

        fglobal = min(Gbuck, Gstress)
        print("Global safety factor: {0:.6g}".format(fglobal))


    def result(self, name):
        """
        Computes a scalar called 'name' that depends
        on the state of the frame. It is used by the plotting functions.
        """
        E = self.theType.YoungModulus
        sigmae = self.theType.sigmae
        A = self.theType.area
        W = self.theType.modulus
        II = self.theType.inertia
        L = self.L
        node0 = self.theNodes[0]
        node1 = self.theNodes[1]

        U = np.array([node0.U[0], node0.U[1], node0.U[2],
                      node1.U[0], node1.U[1], node1.U[2]])

        # initialize values as infinite so as to be updated always
        # any other values might not be modified at the points
        sigmamax = float('-inf')
        sigmamin = float('inf')
        Gstress = float('inf')
        Gbuck = float('inf')

        # calculate the maximum and mininum values of sigma in
        # all points, in the three sections: left, right, and center
        quadPoints = GaussLobattoPoints(3)
        for qp in quadPoints:
            x = 0.5*(qp.x + 1)*L
            B = self.Bmatrix(x)

            [epsilon, kappa] = np.matmul(B, U)
            N = epsilon * E * A
            M = kappa * E * II

            Gstress_k, Gbuck_k = self.safetyFactors(N, M)
            Gstress = min(Gstress, Gstress_k)
            Gbuck = min(Gbuck, Gbuck_k)

            # maximum stress. This needs to be programmed
            sigmamin_k, sigmamax_k = self.stressMinMax(N, M)
            sigmamax = max(sigmamax, sigmamax_k)
            sigmamin = min(sigmamin, sigmamin_k)

        # buckling safefy factor for the whole beam
        print("\nBuckling safety factor: {0:.6g}".format(Gbuck))

        # stress safefy factor for the whole beam
        print("Stress safety factor: {0:.6g}".format(Gstress))

        fglobal = min(Gbuck, Gstress)
        print("Global safety factor: {0:.6g}".format(fglobal))

        # compute some values at the center, for plotting
        qp = quadPoints[1]
        x = 0.5*(qp.x + 1)*L
        B = self.Bmatrix(x)

        [epsilon, kappa] = np.matmul(B, U)
        N = epsilon * E * A
        M = kappa * E * II

        if (name == "epsilon"):
            r = epsilon

        elif (name == "kappa"):
            r = kappa

        elif (name == "N"):
            r = N

        elif (name == "M"):
            r = M

        elif (name == "safety"):
            r = fglobal

        else:
            r = 0.0

        return r


    def rotation(self):
        alpha = self.alpha
        c = math.cos(alpha)
        s = math.sin(alpha)
        r = np.array([[ c, s, 0, 0, 0, 0],
                      [-s, c, 0, 0, 0, 0],
                      [ 0, 0, 1, 0, 0, 0],
                      [ 0, 0, 0, c, s, 0],
                      [ 0, 0, 0,-s, c, 0],
                      [ 0, 0, 0, 0, 0, 1]])
        return r


    def safetyFactors(self, N, M):
        """
        Computes the stress and buckling safety factors.

        Parameters
        ==========
        N: axial force
        M: bending moment

        Output
        ======
        Gstress: stress safety factor. Safe is > 1.
        Gbuck: buckling safety factor. Safe is >1.
        """

        A = self.theType.area
        II = self.theType.inertia
        E = self.theType.YoungModulus
        sigmae = self.theType.sigmae

        Gstress = float('inf')
        Gbuck = float('inf')

        return Gstress, Gbuck


    def stressMinMax(self, N, M):
        """
        Computes the minimum and maximum stresses in a cross section
        of a beam, given the axial force and bending moment it is
        subject to. Note that "min" and "max" DO NOT refer to absolute
        values.

        Parameters
        ==========
        N: axial force
        M: bending moment

        Output
        ======
        sigmamin: The minimum normal stress in the cross section.
        sigmamax: The maximum normal stress in the cross section.
        """
        A = self.theType.area
        W = self.theType.modulus

        sigmamax = 0.0
        sigmamin = 0.0

        return sigmamin, sigmamax


def printIfPositive(title, var):
    print(title, end="")
    if var > 0:
        print(var)
    else:
        print("---")
    return
