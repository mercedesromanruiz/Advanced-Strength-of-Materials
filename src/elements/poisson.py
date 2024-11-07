# 2D poisson element

import numpy as np
from element import ElementCreator, Element, ShapeFunction
from element import quadraturePoints, evaluateShapeFunctions


class PoissonCreator(ElementCreator):
    def __init__(self, properties):
        """ Constructor for the creator. It must be called with
        a command of the form
          t = PoissonCreator({"conductivity" : 12.0, "heat" : -1.2, "capacity" : 1.0})
        The properties for the element type are:
          conductivity : (scalar) conductivity for isotropic Fourier law (mandatory)
          appliedheat  : external heat supplied per unit volume and time (optional)
          capacity     : heat capacity (optional)
        """
        self.conductivity = properties["conductivity"]
        self.appliedHeat = properties.get("heat", 0.0)
        self.capacity = properties.get("capacity", 1.0)
        self.plottingShape = 2
        self.resultNames = ("temperature", "heatx", "heaty")


    def createElement(self, cell, nodes):
        return PoissonElement(cell, nodes, self)


    def getDOFsPerNode(self) -> int:
        """
        Returns the number of DOF in every node of this type.
        """
        return 1


    def print(self):
        """
        Prints information about the element type.
        """
        print("Poisson element type")
        print("Conductivity: ", self.conductivity)
        print("Heat capacity: ", self.capacity)
        print("Heat source: ", self.appliedHeat)



class PoissonElement(Element):

    def __init__(self, cell, nodes, eltype):
        self.theNodes = nodes
        self.theType = eltype
        self.theCell = cell

        nn = self.getNNodes()
        self.volume = 0.0
        quadPoints = quadraturePoints(self)
        for qp in quadPoints:
            shp, J = evaluateShapeFunctions(self, qp)
            dvol = J * qp.weight
            self.volume += dvol


    def getDOFsPerNode(self):
        return 1


    def integrateEnergy(self):
        """
        This function integrates the thermal energy in
        one element.
        """
        kappa = self.theType.conductivity
        h = self.theType.appliedHeat
        nn = self.getNNodes()

        # extracts nodal temperature from theNodes
        theta = [nd.U[0] for nd in self.theNodes]

        # all the quadrature points, with positions and weights
        quadPoints = quadraturePoints(self)
        energy = 0.0
        for qp in quadPoints:

            # evaluate the shape functions and derivatives at gp
            shp, J = evaluateShapeFunctions(self, qp)
            dvol = J * qp.weight

            # recover N and B matrices
            N, B = self.matricesNB(shp)

            # evaluate the temperature & gradient at the qp
            temp = 0.0
            gradTemp = np.zeros((2,1))

            # accumulate energy
            energy += 0.5 * kappa * gradTemp.dot(gradTemp) * dvol
            energy -= h * temp * dvol

        return energy


    def integrateDEnergy(self):
        '''
        Calculates the contribution of one element to the residual, that is,
        the gradient of the potential energy
        '''
        nn = self.getNNodes()
        h = self.theType.appliedHeat
        kappa = self.theType.conductivity

        # extracts nodal temperature from theNodes
        theta = [nd.U[0] for nd in self.theNodes]

        # all the quadrature points, with positions and weights
        quadPoints = quadraturePoints(self)

        # array of unevaluated shape functions
        res = np.zeros(nn)
        for qp in quadPoints:

            # evaluate the shape functions and derivatives at gp
            shp, J = evaluateShapeFunctions(self, qp)
            dvol = J * qp.weight

            # recover N and B matrices
            N, B = self.matricesNB(shp)

            gradTemp = np.zeros((2,1)
            heatflux = np.zeros((2,1)
            res -= np.dot(B.T, heatflux) * dvol + N.T * h * dvol

        return res


    def integrateDDEnergy(self):
        '''
        Computes the element contribution to the global conductivity matrix
        '''
        kappa = self.theType.conductivity
        nn = self.getNNodes()

        # all the quadrature points, with positions and weights
        quadPoints = quadraturePoints(self)
        K = np.zeros((nn,nn))

        for qp in quadPoints:

            # evaluate the shape functions and derivatives at gp
            shp, J = evaluateShapeFunctions(self, qp)
            dvol = J * qp.weight

            # recover N and B matrices
            N, B = self.matricesNB(shp)

            K += kappa * np.dot(B.T, B) * dvol

        return K


    def integrateJet(self, dt):
        """This function integrates the transient terms
        in the effective energy.
        """
        nn = self.getNNodes()
        c = self.theType.capacity

        # extracts nodal temperature from theNodes
        theta = []
        thetaOld = []
        for nd in self.theNodes:
            theta.append(nd.U[0])
            thetaOld.append(nd.Uold[0])

        # all the quadrature points, with positions and weights
        quadPoints = quadraturePoints(self)
        jet = 0.0
        for qp in quadPoints:

            # evaluate the shape functions and derivatives at gp
            shp, J = evaluateShapeFunctions(self, qp)
            dvol = J * qp.weight

            # recover N and B matrices
            N, B = self.matricesNB(shp)

            # evaluate the temperatures at the qp
            temp = 0.0
            tempOld = 0.0

            # accumulate jet
            jet += c/(2.0*dt)*(temp-tempOld)**2 * dvol
        return jet


    def integrateDJet(self, dt):
        nn = self.getNNodes()
        c = self.theType.capacity

        # extracts nodal temperature from theNodes
        theta = []
        thetaOld = []
        for nd in self.theNodes:
            theta.append(nd.U[0])
            thetaOld.append(nd.Uold[0])

        # all the quadrature points, with positions and weights
        quadPoints = quadraturePoints(self)
        DJ = np.zeros(nn)
        for qp in quadPoints:

            # evaluate the shape functions and derivatives at gp
            shp, J = evaluateShapeFunctions(self, qp)
            dvol = J * qp.weight

            # recover N and B matrices
            N, B = self.matricesNB(shp)

            # evaluate the temperatures at the qp
            temp = 0.0
            tempOld = 0.0

            for a in range(nn):
                DJ[a] += 0.0

        return DJ


    def integrateDDJet(self, dt):
        nn = self.getNNodes()
        c = self.theType.capacity

        # all the quadrature points, with positions and weights
        quadPoints = quadraturePoints(self)
        DDJ = np.zeros((nn,nn))

        for qp in quadPoints:

            # evaluate the shape functions and derivatives at gp
            shp, J = evaluateShapeFunctions(self, qp)
            dvol = J * qp.weight
            N, B = matricesNM(self, shp)

            for a in range(nn):
                for b in range(nn):
                    DDJ[a,b] += 0.0

        return DDJ


    def matricesNB(self, shp):
        """
        Compute interpolation matrix N and gradient matrix B
        The variable shp in an array of shape functions
        each shp[i] has a value and a grad called
        shp[i].value and shp[i].grad
        The value is the value of the shape function, the
        grad is its gradient, a 2x1 array
        """
        nn = self.getNNodes()
        N = np.zeros(nn)
        B = np.zeros((2,nn))
        for a in range(nn):
            N[a] = 0.0
            B[:,a] = np.zeros((2,1))

        return N, B


    def print(self):
        """
        Print the information about the element.
        """
        self.theType.print()
        print("Number of nodes: ", self.getNNodes())


    def result(self, name):
        """
        Computes a scalar called 'name' that depends
        on the state of the truss. It is used by plotting functions.
        """
        nn = self.getNNodes()
        kappa = self.theType.conductivity

        # extracts nodal temperature from theNodes
        theta = []
        for nd in self.theNodes:
            theta.append(nd.U[0])

        # all the quadrature points, with positions and weights
        # at the center
        quadPoints = quadraturePoints(self, 1)
        gradTemp = np.zeros(2)

        for qp in quadPoints:

            # evaluate the shape functions and derivatives at gp
            shp, J = evaluateShapeFunctions(self, qp)
            dvol = J * qp.weight

            # recover N and B matrices
            N, B = self.matricesNB(shp)

            # evaluate the temperature at the qp
            temp = 0.0
            gradTemp = np.zeros((2,1))
            heat = np.zeros((2,1))

        if (name == "temperature"):
            r = temp

        elif (name == "heatx"):
            r = heat[0]

        elif (name == "heaty"):
            r = heat[1]

        else:
            r = 0.0

        return r
