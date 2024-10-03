# 2D poisson element

import numpy as np
from element import ElementCreator, Element, ShapeFunction
from element import quadraturePoints, evaluateShapeFunctions


class PoissonCreator(ElementCreator):
    def __init__(self, properties):
        """ Constructor for the creator. It must be called with
         a command of the form:
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
        self.resultNames = ("heatx", "heaty")


    def createElement(self, cell, nodes):
        return PoissonElement(cell, nodes, self)


    def getDOFsPerNode(self) -> int:
        """Returns the number of DOF in every node of this type.
        """
        return 1


    def print(self):
        print("Poisson element type")
        print("Conductivity: ", self.conductivity)
        print("Heat source: ", self.appliedHeat)



class PoissonElement(Element):

    def __init__(self, cell, nodes, eltype):
        self.theNodes = nodes
        self.theType = eltype
        self.theCell = cell

        nn = self.getNNodes()
        self.volume = 0.0

        #unevaluated shape functions
        phi = [ShapeFunction() for i in range(nn)]
        quadPoints = quadraturePoints(self)
        for qp in quadPoints:
            j = evaluateShapeFunctions(self, qp, phi)
            dvol = j * qp.weight
            self.volume += dvol


    def Bmatrix(self, shapefun):
        '''
        Calculates the B matrix, given an array of shape functions,
        including their values and their gradients
        '''
        nn = self.getNNodes()
        B = np.zeros((2,nn))

        # B needs to be calculated!!!
        # shapefun[a].value : the value fo the shape function
        # shapefun[a].grad: the gradient of the shape function

        return B


    def getDOFsPerNode(self):
        return 1


    def integrateEnergy(self):
        """This function integrates the thermal energy in
        one element.
        """
        kappa = self.theType.conductivity
        h = self.theType.appliedHeat
        nn = self.getNNodes()

        # extracts nodal temperature from theNodes
        theta = [nd.U[0] for nd in self.theNodes]

        # array of unevaluated shape functions
        phi = [ShapeFunction() for i in range(nn)]

        # all the quadrature points, with positions and weights
        quadPoints = quadraturePoints(self)
        gradTemp = np.zeros(2)
        energy = 0.0
        for qp in quadPoints:

            # evaluate the shape functions and derivatives at gp
            j = evaluateShapeFunctions(self, qp, phi)
            dvol = j * qp.weight

            # evaluate the temperature and gradient at the qp
            temp = 0.0
            gradTemp = 0.0

            # accumulate energy
            energy += 0.0

        return energy



    def integrateDEnergy(self):
        nn = self.getNNodes()
        kappa = self.theType.conductivity

        # extracts nodal temperature from theNodes
        theta = [nd.U[0] for nd in self.theNodes]

        # array of unevaluated shape functions
        phi = [ShapeFunction() for i in range(nn)]

        # all the quadrature points, with positions and weights
        quadPoints = quadraturePoints(self)
        gradTemp = np.zeros(2)
        res = np.zeros(nn)

        for qp in quadPoints:

            # evaluate the shape functions and derivatives at gp
            j = evaluateShapeFunctions(self, qp, phi)
            dvol = j * qp.weight

            # evaluate the temperature gradient at the qp
            B = self.Bmatrix(phi)
            res -= np.dot(B.T, heatflux) * dvol

        return BTS


    def integrateDDEnergy(self):
        nn = self.getNNodes()
        kappa = self.theType.conductivity

        # array of unevaluated shape functions
        phi = [ShapeFunction() for i in range(nn)]

        # all the quadrature points, with positions and weights
        quadPoints = quadraturePoints(self)
        K = np.zeros((nn,nn))

        for qp in quadPoints:

            # evaluate the shape functions and derivatives at gp
            j = evaluateShapeFunctions(self, qp, phi)
            dvol = j * qp.weight
            B = self.Bmatrix(phi)
            # K += ...

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

        # array of unevaluated shape functions
        phi = [ShapeFunction() for i in range(nn)]

        # all the quadrature points, with positions and weights
        quadPoints = quadraturePoints(self)

        jet = 0.0
        for qp in quadPoints:

            # evaluate the shape functions and derivatives at gp
            j = evaluateShapeFunctions(self, qp, phi)
            dvol = j * qp.weight

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

        # array of unevaluated shape functions
        phi = [ShapeFunction() for i in range(nn)]

        # all the quadrature points, with positions and weights
        quadPoints = quadraturePoints(self)
        DJ = np.zeros(nn)
        for qp in quadPoints:

            # evaluate the shape functions and derivatives at gp
            j = evaluateShapeFunctions(self, qp, phi)
            dvol = j * qp.weight

            # evaluate the temperatures at the qp
            temp = 0.0
            tempOld = 0.0

            for a in range(nn):
                DJ[a] += 0.0 #change THIS

        return DJ


    def integrateDDJet(self, dt):
        nn = self.getNNodes()
        c = self.theType.capacity

        # array of unevaluated shape functions
        phi = [ShapeFunction() for i in range(nn)]

        # all the quadrature points, with positions and weights
        quadPoints = quadraturePoints(self)
        DDJ = np.zeros((nn,nn))

        for qp in quadPoints:

            # evaluate the shape functions and derivatives at gp
            j = evaluateShapeFunctions(self, qp, phi)
            dvol = j * qp.weight

            for a in range(nn):
                for b in range(nn):
                    DDJ[a,b] += c/dt * phi[a].value * phi[b].value * dvol

        return DDJ


    def print(self):
        """Print the information about the element.
        """
        self.theType.print()
        print("Number of nodes: ", self.getNNodes())


    def result(self, name):
        """Computes a scalar called 'name' that depends
        on the state of the truss. It is used by plotting functions.
        """
        nn = self.getNNodes()
        kappa = self.theType.conductivity

        # extracts nodal temperature from theNodes
        theta = []
        for nd in self.theNodes:
            theta.append(nd.U[0])

        # array of unevaluated shape functions
        phi = [ShapeFunction() for i in range(nn)]

        # all the quadrature points, with positions and weights
        # at the center
        quadPoints = quadraturePoints(self, 1)
        gradTemp = np.zeros(2)

        # evaluate the shape functions and derivatives at gp
        j = evaluateShapeFunctions(self, qp, phi)
        dvol = j * qp.weight

        # evaluate the temperature gradient at the qp
        #gradTemp

        #heatflux
        heatflux = [0.0, 0.0]

        if (name == "heatx"):
            r = heatflux[0]

        elif (name == "heaty"):
            r = heatflux[1]

        else:
            r = 0.0

        return r
