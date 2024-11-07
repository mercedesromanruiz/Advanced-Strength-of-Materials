# 2d/3d mechanical problems with elastic materials

import numpy as np
import numpy.linalg as LA

from element import ElementCreator, Element, ShapeFunction
from element import quadraturePoints, evaluateShapeFunctions


class MechanicalCreator(ElementCreator):
    def __init__(self, properties):
        """
        Set up an element creator, initializing
        the data for the elastic material and the volumetric
        loading.
           young: Young's modulus
           poisson: Poisson's ratio
           bx : x-component of the volumetric load
           by : y-component of the volumetric load
           bz : z-component of the volumetric load
        """
        self.young = properties.get("young", 1.0)
        self.poisson = properties.get("poisson", 0.0)
        self.b = np.zeros(3)
        self.b[0] = properties.get("bx", 0.0)
        self.b[1] = properties.get("by", 0.0)
        self.b[2] = properties.get("bz", 0.0)

        self.plottingShape = 3
        self.resultNames = ("vonmises", "tresca",
                            "pressure", "volstrain", "max strain", "min strain",
                            "sigmaxx", "sigmayy", "sigmazz", "sigmayz", "sigmaxz",
                            "sigmaxy")


    def createElement(self, cell, nodes):
        self.plottingShape = len(nodes[0].coordinates())
        return MechanicalElement(cell, nodes, self)


    def getDOFsPerNode(self) -> int:
        """Returns the number of DOF in every node of this type.
        """
        dim = self.plottingShape
        return dim


    def print(self):
        dim = self.getDOFsPerNode()
        print("Mechanical element type")
        print("Young modulus:   ", self.young)
        print("Poisson's ratio: ", self.poisson)
        print("Body force : ", [self.b[i] for i in range(dim)])


class MechanicalElement(Element):

    def __init__(self, cell, nodes, eltype):
        self.theNodes = nodes
        self.theType = eltype
        self.theCell = cell
        self.volume = 0.0

    def buildDisplacementVector(self):
        dim = self.theType.getDOFsPerNode()
        nen = self.getNNodes()
        disp = np.zeros((nen,dim))
        for a, nd in enumerate(self.theNodes):
            if (dim == 2):
                disp[a] = [ nd.U[0], nd.U[1] ]
            else:
                disp[a] = [ nd.U[0], nd.U[1], nd.U[2] ]
        return disp


    def computeUhGradUh(self, phi, disp):
        """
        Calculate the displacement u at the gauss points and its gradient
        Input:
         phi: an array of shape functions and their gradients
         disp: an np vector of nodal displacments [U1x, U1y, U2x, U2y, ...]
        Output:
         uh(x)
         Gradient[uh(x)]
        """
        dim = self.theType.getDOFsPerNode()
        nen = len(phi)
        Uh = np.zeros(dim)
        gradUh = np.zeros((dim, dim))

        for a in range(nen):
            Uh     += phi[a].value * disp[a]
            gradUh += np.outer(disp[a], phi[a].grad)

        return Uh, gradUh


    def getDOFsPerNode(self):
        return self.theType.getDOFsPerNode()


    def integrateEnergy(self):
        """
        This function integrates the mechanical energy in
        one element. This includes the stored energy from the
        deformation and the potential energy of the external body
        forces (bx,by,bz).
        """
        dim = self.theType.getDOFsPerNode()
        nen = self.getNNodes()
        E  = self.theType.young
        nu = self.theType.poisson
        ll = nu*E/(1.0-2.0*nu)/(1.0+nu)
        mu = E/2.0/(1.0+nu)
        bf = np.array( [self.theType.b[i] for i in range(dim)])

        # extracts nodal displacements from theNodes
        disp = self.buildDisplacementVector()

        # all the quadrature points, with positions and weights
        quadPoints = quadraturePoints(self)

        strain = np.zeros((dim, dim))
        energy = 0.0
        for qp in quadPoints:

            # evaluate the shape functions and derivatives at gp
            phi, J = evaluateShapeFunctions(self, qp)
            dvol = J * qp.weight

            # evaluate the displacement u and grad(u) at the qp
            Uh, gradUh = self.computeUhGradUh(phi, disp)

            # evaluate strain and its trace
            strain = 0.5 * (gradUh + gradUh.T)
            traceE = 0.0
            for i in range(dim): traceE += strain[i,i]

            # accumulate internal and external energy
            energy += ??
        return energy


    def integrateDEnergy(self):
        dim = self.theType.getDOFsPerNode()
        nen = self.getNNodes()
        E  = self.theType.young
        nu = self.theType.poisson
        ll = nu*E/(1.0-2.0*nu)/(1.0+nu)
        mu = E/2.0/(1.0+nu)
        bf = np.array( [self.theType.b[i] for i in range(dim)])

        # extracts nodal displacements from the nodes
        disp = self.buildDisplacementVector()

        # all the quadrature points, with positions and weights
        quadPoints = quadraturePoints(self)

        strain = np.zeros((dim,dim))
        stress = np.zeros((dim,dim))
        eye = np.identity(dim)
        F = np.zeros(nen*dim)
        for qp in quadPoints:

            # evaluate the shape functions and derivatives at gp
            phi, J = evaluateShapeFunctions(self, qp)
            dvol = J * qp.weight

            # evaluate the displacement u and grad(u) at the qp
            Uh, gradUh = self.computeUhGradUh(phi, disp)

            # evaluate strain and its trace
            strain = 0.5 * (gradUh + gradUh.T)
            traceE = 0.0
            for i in range(dim): traceE += strain[i,i]

            # evaluate stress tensor
            stress = ll * traceE * eye + 2.0 * mu * strain

            # accumulate internal forces
            for a in range(nen):
                a0 = dim*a
                a1 = a0 + dim
                F[a0:a1] += ??

        return F


    def integrateDDEnergy(self):
        dim = self.theType.getDOFsPerNode()
        nen = self.getNNodes()
        E  = self.theType.young
        nu = self.theType.poisson
        ll = nu*E/(1.0-2.0*nu)/(1.0+nu)
        mu = E/2.0/(1.0+nu)
        eye = np.identity(dim)
        K = np.zeros((nen*dim, nen*dim))

        # all the quadrature points, with positions and weights
        quadPoints = quadraturePoints(self)

        for qp in quadPoints:

            # evaluate the shape functions and derivatives at gp
            phi, J = evaluateShapeFunctions(self, qp)
            dvol = J * qp.weight

            for a in range(nen):
                a0 = dim*a
                a1 = a0 + dim
                for b in range(nen):
                    b0 = dim*b
                    b1 = b0 + dim

                    K[a0:a1,b0:b1] += ll * np.outer(phi[a].grad, phi[b].grad) * dvol
                    K[a0:a1,b0:b1] += mu * np.dot(phi[a].grad, phi[b].grad) * eye * dvol
                    K[a0:a1,b0:b1] += mu * np.outer(phi[b].grad, phi[a].grad) * dvol

        return K


    def print(self):
        """Print the information about the element.
        """
        self.theType.print()
        print("Number of nodes: ", self.getNNodes())


    def result(self, name):
        """Computes a scalar called 'name' that depends
        on the state of the solid. It is used by plotting functions.
        """
        dim = self.theType.getDOFsPerNode()
        nen = self.getNNodes()
        E  = self.theType.young
        nu = self.theType.poisson
        ll = nu*E/(1.0-2.0*nu)/(1.0+nu)
        mu = E/2.0/(1.0+nu)
        eye = np.identity(dim)

        # extracts nodal displacements from theNodes
        disp = []
        for a,nd in enumerate(self.theNodes):
            if (dim == 2):
                disp.append(np.array([ nd.U[0], nd.U[1] ]))
            else:
                disp.append(np.array([ nd.U[0], nd.U[1], nd.U[2] ]))

        # all the quadrature points, with positions and weights
        quadPoints = quadraturePoints(self, 1)
        qp = quadPoints[0]

        gradU = np.zeros((dim,dim))
        strain = np.zeros((dim,dim))
        stress = np.zeros((dim,dim))
        u = np.zeros(dim)

        # evaluate the shape functions and derivatives at gp
        phi, J = evaluateShapeFunctions(self, qp)
        dvol = J * qp.weight

        # evaluate the displacement at the qp
        u.fill(0.0)
        gradU.fill(0.0)
        for a in range(nen):
            u += phi[a].value * disp[a]
            gradU += np.dot(phi[a].grad, disp[a].T)

        # evaluate strain and stress
        strain = (gradU + gradU.T) * 0.5
        traceE = 0.0
        for i in range(dim):
            traceE += strain[i,i]

        stress = ll * traceE * eye + 2.0*mu*strain

        if (name == "vonmises"):
            r = 0.0

        elif (name == "tresca"):
            r = 0.0

        elif (name == "pressure"):
            p = 0.0
            for i in range(dim):
                p -= stress[i,i]
            r = p/dim

        elif (name == "volstrain"):
            theta = 0.0
            for i in range(dim):
                theta += strain[i,i]
            r = theta

        elif (name == "maxstrain"):
            ev = LA.eigvals(strain)
            evs = ev.argsort()[::-1]
            r = evs[dim]

        elif (name == "minstrain"):
            ev = LA.eigvals(strain)
            evs = ev.argsort()[::-1]
            r = evs[0]

        elif (name == "sigmaxx"):
            r = stress[0,0]

        elif (name == "sigmayy"):
            r = stress[1,1]

        elif (dim > 2 and name == "sigmazz"):
            r = stress[2,2]

        elif (dim > 2 and name == "sigmayz"):
            r = stress[1,2]

        elif (dim > 2 and name == "sigmaxz"):
            r = stress[0,2]

        elif (name == "sigmaxy"):
            r = stress[0,1]

        else:
            r = 0.0

        return r
