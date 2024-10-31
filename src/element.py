import numpy as np
import math


class ElementCreator:
    """This class embodies the producer pattern. Each ElementCreator
    is later used to spawn many elements of the same type, with
    the same parameters.
    """

    def __init__(self):
        self.plottingShape = 0
        self.resultNames = {}
        pass


    def createElement(self):
        """This is a pure virtual class, the main task of every
        ElementCreator.
        """
        pass


    def getDOFsPerNode(self) -> int:
        """This is a pure virtual class. It should never be called
        """
        print("Error calling getDOFsPerNode in parent class.")
        return 0



class Element:
    """
    This is the parent class for all elements in a finite element
    models.
    """

    def __init__(self, cell, eltype):
        self.theCell = cell
        self.theType = eltype
        self.volume = 0.0


    def getDOFs(self):
        dofs = []
        dpn = self.theType.getDOFsPerNode()
        for n in self.theNodes:
            for i in range(dpn):
                dofs.append(n.DOFS[i])

        return dofs


    def getNDOFs(self):
        return self.getNNodes() * self.theType.getDOFsPerNode()


    def getNNodes(self):
        return self.theCell.getNVertices()


    def integrateBmatrix(self):
        print("IntegrateBmatrix should not be used")
        pass


    def integrateEnergy(self) -> float:
        en = 0.0
        return en


    def integrateDEnergy(self):
        print("IntegrateDEnergy parent class")
        pass


    def integrateDDEnergy(self):
        pass


    def integrateJet(self, dt) -> float:
        jet = 0.0
        return jet


    def integrateDJet(self, dt) -> float:
        return 0.0


    def integrateDDJet(self, dt) -> float:
        return 0.0


    def integrateVolume(self):
        return 0.0


    def numberDOFs(self, nextdof):
        for n in self.theNodes:
            nextdof = n.numberConsecutiveDOFs(nextdof)
        return nextdof


    def print(self):
        print("in parent print element")
        pass


    def result(self, name):
        print("in parent element, function result")
        pass



class ShapeFunction:
    """A shape function is an objects that holds a value
    and the derivatives with respect to three coordinates (the
    grad-ient). Maybe one or more of the latter can be empty.
    """
    def __init__(self):
        self.value = 0.0
        self.grad = np.array(0)


def evaluateShapeFunctions(el, qp):
    """
    Evaluates a list of shape functions for an element and
    a quadrature point, and also computes the jacobian j of
    the element transformation.
    """
    nn = el.getNNodes()

    shp = [ShapeFunction() for i in range(nn)]

    for a in range(nn):
        shp[a].grad = np.zeros(el.theType.plottingShape)

    if (el.theType.plottingShape == 1):
        j = evaluateShapeFunctions1d(shp, el.theNodes, qp)
    elif (el.theType.plottingShape == 2):
        j = evaluateShapeFunctions2d(shp, el.theNodes, qp)
    elif (el.theType.plottingShape == 3):
        j = evaluateShapeFunctions3d(shp, el.theNodes, qp)

    return shp, j


def evaluateShapeFunctions1d(theShp, theNodes, qp):
    """ Evaluates a list of shape functions for 1-dimensional
    problems.
    """
    nn = len(theShp)
    xi = qp.x

    if (nn == 2):
        nodei = theNodes[0]
        nodej = theNodes[1]

        r = nodej.coordinates() - nodei.coordinates()
        L = math.sqrt(np.dot(r,r))
        iL = 1.0/L
        j = 0.5*L

        theShp[0].value = 0.5 * (1.0 - xi)
        theShp[1].value = 0.5 * (1.0 + xi)

        theShp[0].grad[0] = -iL
        theShp[1].grad[0] = +iL


    # quadratic lagrange polynomials. Numbering 1 --- 3 --- 2
    elif (nn == 3):

        theShp[0].value = 0.5 * xi * (xi - 1.0);
        theShp[1].value = 0.5 * xi * (xi + 1.0);
        theShp[2].value = 1.0 - xi * xi;

        # natural derivatives
        theShp[0].grad[0] = xi - 0.5;
        theShp[1].grad[0] = xi + 0.5;
        theShp[2].grad[0] = -2.0*xi;

        # jacobian
        v = np.zeros(3)
        for a in range(nn):
            v += theShp[a].grad[0] * theNodes[a].coordinates();
        j = math.sqrt(v.dot(v))

        for a in range(nn):
           theShp[a].grad[0] /= j;

    else:
        print("1d shape function not programmed for this n nodes")

    return j



def evaluateShapeFunctions2d(theShp, theNodes, qp):

    nn = len(theNodes)
    xi  = qp.x
    eta = qp.y
    dN = []

    if (nn == 3):
        zeta = 1.0 - xi - eta

        theShp[0].value = xi
        theShp[1].value = eta
        theShp[2].value = zeta

        dN.append( np.array([ 1.0,  0.0]) )
        dN.append( np.array([ 0.0,  1.0]) )
        dN.append( np.array([-1.0, -1.0]) )

        j = pushShapeFunctions(dN, theShp, theNodes)

    elif (nn == 4):

        theShp[0].value = 0.25 * (1.0 + xi) * (1.0 + eta)
        theShp[1].value = 0.25 * (1.0 - xi) * (1.0 + eta)
        theShp[2].value = 0.25 * (1.0 - xi) * (1.0 - eta)
        theShp[3].value = 0.25 * (1.0 + xi) * (1.0 - eta)

        dN.append( np.array([ 0.25*(1.0+eta),  0.25*(1.0+xi)]) )
        dN.append( np.array([-0.25*(1.0+eta),  0.25*(1.0-xi)]) )
        dN.append( np.array([-0.25*(1.0-eta), -0.25*(1.0-xi)]) )
        dN.append( np.array([ 0.25*(1.0-eta), -0.25*(1.0+xi)]) )

        j = pushShapeFunctions(dN, theShp, theNodes)

    else:
        print("2d shape function not programmed for this n nodes")

    return j


def evaluateShapeFunctions3d(theShp, theNodes, qp):
    """ evaluate the shape functions of a 3d element.
    Only programmed for tetrahedra yet.
    """
    nn = len(theNodes)
    xi = np.zeros(4)

    if (nn == 4):

        x = np.zeros(4)
        y = np.zeros(4)
        z = np.zeros(4)
        for i, nd in enumerate(theNodes):
            x[i] = nd.coordinates()[0]
            y[i] = nd.coordinates()[1]
            z[i] = nd.coordinates()[2]

        perm = (0, 1, 2, 3, 0, 1, 2)
        a = np.zeros((4,4))

        for i in range(4):

            j = perm[i+1]
            k = perm[i+2]
            l = perm[i+3]

            a[i,1] = y[k]*z[l] - y[l]*z[k] \
            +        y[l]*z[j] - y[j]*z[l] \
            +        y[j]*z[k] - y[k]*z[j]

            a[i,2] = x[k]*z[l] - x[l]*z[k] \
            +        x[l]*z[j] - x[j]*z[l] \
            +        x[j]*z[k] - x[k]*z[j]

            a[i,3] = x[k]*y[l] - x[l]*y[k] \
            +        x[l]*y[j] - x[j]*y[l] \
            +        x[j]*y[k] - x[k]*y[j]

            a[i,j] = -a[i,j];
            a[i,l] = -a[i,l];

        det = x[0]*a[0,1] + x[1]*a[1,1] + x[2]*a[2,1] + x[3]*a[3,1];
        xi[0] = qp.x
        xi[1] = qp.y
        xi[2] = qp.z
        xi[3] = 1.0 - xi[0] - xi[1] - xi[2];

        idet = 1.0/det;
        for i in range(4):
            theShp[i].value = xi[i]
            theShp[i].grad[0] = a[i,1]*idet
            theShp[i].grad[1] = a[i,2]*idet
            theShp[i].grad[2] = a[i,3]*idet

        j = det/6.0;

    elif (nn == 8):

        xisgn = (+1.0, -1.0, -1.0, +1.0, +1.0, -1.0, -1.0, +1.0)
        etsgn = (+1.0, +1.0, -1.0, -1.0, +1.0, +1.0, -1.0, -1.0)
        zesgn = (-1.0, -1.0, -1.0, -1.0, +1.0, +1.0, +1.0, +1.0)

        # reference element coordinates
        xi  = qp.x
        eta = qp.y
        zeta = qp.z

        # the value of the shape functions at the quad point
        for a in range(8):
            theShp[a].value = 0.125*(1.0+xisgn[a]*xi)*(1.0+etsgn[a]*eta)*(1.0+zesgn[a]*zeta)


        # the shape function natural derivatives
        Nd = np.zeros([8,3])
        for a in range(8):
            Nd[a,0] = 0.125 *        xisgn[a]      * (1.0 + etsgn[a]*eta)  * (1.0 + zesgn[a]*zeta)
            Nd[a,1] = 0.125 * (1.0 + xisgn[a]*xi ) *        etsgn[a]       * (1.0 + zesgn[a]*zeta)
            Nd[a,2] = 0.125 * (1.0 + xisgn[a]*xi ) * (1.0 + etsgn[a]*eta ) *        zesgn[a]

        j = pushShapeFunctions(Nd, theShp, theNodes)


    else:
        print("3d shape function not programmed for this n nodes", nn)

    return j


def pushShapeFunctions(dN, theShp, theNodes):
    """Pushes the isoparametric gradients to the current gradients
    and stores the latter in the shape function"""

    dim = len(theNodes[0].coordinates())
    J = np.zeros((dim,dim))
    nn = len(theNodes)
    for a in range(nn):
        J = J + np.outer(theNodes[a].coordinates(), dN[a])

    j = np.linalg.det(J)
    Jit = np.linalg.inv(J).T

    for a in range(nn):
        theShp[a].grad = Jit.dot(dN[a])

    return j



class Quadpoint:

    def __init__(self, x=0.0, y=0.0, z=0.0, w=0.0):
        self.weight = w
        self.x = x
        self.y = y
        self.z = z

    def print(self):
        print("Quadradure point")
        print("Weight: ", self.weight)
        print("Coordinates: ", self.x, self.y, self.z)


def quadraturePoints(e, nqp=0):
    shape = e.theType.plottingShape

    if shape == 1:
        qps = quadraturePoints1D(nqp)

    elif shape == 2 and e.getNNodes() == 3:
        qps = quadraturePointsTriangle(nqp)

    elif shape == 2 and e.getNNodes() == 4:
        qps = quadraturePointsQuad(nqp)

    elif shape == 3 and e.getNNodes() == 4:
        qps = quadraturePointsTetrahedron(nqp)

    elif shape == 3 and e.getNNodes() == 8:
        qps = quadraturePointsHexahedron(nqp)

    else:
        print("Quadrature rule not defined")

    return qps


def quadraturePoints1D(nqp):
    if (nqp == 0):
        nqp = 1
    qps = [Quadpoint() for i in range(nqp)]

    if (nqp == 1):
        qps[0] = Quadpoint(0.0 ,0.0 ,0.0, 2.0)

    elif nqp == 2:
        sq3 = 1.0/math.sqrt(3.0)
        qps[0] = Quadpoint(-sq3, 0.0, 0.0, 1.0)
        qps[1] = Quadpoint(+sq3, 0.0, 0.0, 1.0)

    elif nqp == 3:
        sq35 = math.sqrt(3.0/5.0)
        qps[0] = Quadpoint(-sq35, 0.0, 0.0, 5.0/9.0)
        qps[1] = Quadpoint(  0.0, 0.0, 0.0, 8.0/9.0)
        qps[2] = Quadpoint(+sq35, 0.0, 0.0, 5.0/9.0)

    elif nqp == 4:
        sq1 = math.sqrt(3.0/7.0 - 2.0/7.0*math.sqrt(6.0/5.0))
        w1  = (18.0 + math.sqrt(30.0))/36.0
        sq2 = math.sqrt(3.0/7.0 + 2.0/7.0*math.sqrt(6.0/5.0))
        w2  = (18.0 - math.sqrt(30.0))/36.0

        qps[0] = Quadpoint(-sq1, 0.0, 0.0, w1)
        qps[1] = Quadpoint(+sq1, 0.0, 0.0, w1)
        qps[2] = Quadpoint(-sq2, 0.0, 0.0, w2)
        qps[3] = Quadpoint(+sq2, 0.0, 0.0, w2)

    elif nqp == 5:
        sq1 = 1.0/3.0*math.sqrt(5.0-2.0*math.sqrt(10.0/7.0))
        w1  = (322.0 + 13.0*math.sqrt(70.0))/900.0
        sq2 = 1.0/3.0*math.sqrt(5.0+2.0*math.sqrt(10.0/7.0))
        w2  = (322.0 - 13.0*math.sqrt(70.0))/900.0

        qps[0] = Quadpoint(-sq1, 0.0, 0.0, w1)
        qps[1] = Quadpoint(+sq1, 0.0, 0.0, w1)
        qps[2] = Quadpoint(-sq2, 0.0, 0.0, w2)
        qps[3] = Quadpoint(+sq2, 0.0, 0.0, w2)
        qps[4] = Quadpoint( 0.0, 0.0, 0.0, 128.0/225.0)

    else:
        print("No quadrature rule for shp1 and ngp: ", nqp)

    return qps



def quadraturePointsTriangle(nqp):
    qps = []
    if nqp == 1:
        qps.append( Quadpoint(1.0/3.0, 1.0/3.0, 0.0, 0.5) )

    elif (nqp == 0 or nqp == 3):
        cst = (0.1666666666667, 0.6666666666667)
        cstw = 1.0/6.0

        qps.append( Quadpoint( cst[0], cst[0], 0.0, cstw) )
        qps.append( Quadpoint( cst[1], cst[0], 0.0, cstw) )
        qps.append( Quadpoint( cst[0], cst[1], 0.0, cstw) )

    else:
        print("No quadrature rule for shp1 and ngp: ", nqp)

    return qps


def quadraturePointsQuad(nqp):
    qps = []
    if nqp == 1:
        qps.append( Quadpoint(0.0 ,0.0 ,0.0, 4.0) )

    elif (nqp == 0 or nqp == 4):
        sq3 = 1.0/math.sqrt(3.0)
        qps.append(Quadpoint(+sq3, +sq3, 0.0, 1.0))
        qps.append(Quadpoint(-sq3, +sq3, 0.0, 1.0))
        qps.append(Quadpoint(-sq3, -sq3, 0.0, 1.0))
        qps.append(Quadpoint(+sq3, -sq3, 0.0, 1.0))

    else:
        print("No quadrature rule for Quad and points: ", nqp)

    return qps


def quadraturePointsTetrahedron(nqp):
    qps = []

    if nqp == 1:
        qps.append( Quadpoint(0.25, 0.25, 0.25, 1.0) )

    elif (nqp == 0 or nqp == 4):

        tet4 =(0.5854101966249685, 0.138196601125015, 0.1381966011250105, 0.1381966011250105)
        eps4 = (0, 1, 2, 3, 0, 1)

        for a in range(4):
            qps.append(
                Quadpoint(tet4[eps4[a]], tet4[eps4[a+1]], tet4[eps4[a+2]], 0.25))

    else:
        print("No quadrature rule for Tet and points: ", nqp)

    return qps


def quadraturePointsHexahedron(nqp):
    qps = []
    if nqp == 1:
        qps.append( Quadpoint(0.0 ,0.0 ,0.0, 8.0) )

    elif (nqp == 0 or nqp == 8):
        sq3 = 1.0/math.sqrt(3.0)
        qps.append(Quadpoint(+sq3, +sq3, +sq3, 1.0))
        qps.append(Quadpoint(-sq3, +sq3, +sq3, 1.0))
        qps.append(Quadpoint(-sq3, -sq3, +sq3, 1.0))
        qps.append(Quadpoint(+sq3, -sq3, +sq3, 1.0))
        qps.append(Quadpoint(+sq3, +sq3, -sq3, 1.0))
        qps.append(Quadpoint(-sq3, +sq3, -sq3, 1.0))
        qps.append(Quadpoint(-sq3, -sq3, -sq3, 1.0))
        qps.append(Quadpoint(+sq3, -sq3, -sq3, 1.0))

    else:
        print("No quadrature rule for hex and points: ", nqp)

    return qps


def GaussLobattoPoints(nqp):
    qps = [Quadpoint() for i in range(nqp)]

    if (nqp == 1):
        qps[0] = Quadpoint(0.0 ,0.0 ,0.0, 2.0)

    elif nqp == 2:
        qps[0] = Quadpoint(-1.0, 0.0, 0.0, 1.0)
        qps[1] = Quadpoint(+1.0, 0.0, 0.0, 1.0)

    elif nqp == 3:
        qps[0] = Quadpoint(-1.0, 0.0, 0.0, 1.0/3.0)
        qps[1] = Quadpoint( 0.0, 0.0, 0.0, 4.0/3.0)
        qps[2] = Quadpoint(+1.0, 0.0, 0.0, 1.0/3.0)

    elif nqp == 4:
        sq15 = math.sqrt(1.0/5.0)

        qps[0] = Quadpoint( -1.0, 0.0, 0.0, 1.0/6.0)
        qps[1] = Quadpoint(-sq15, 0.0, 0.0, 5.0/6.0)
        qps[2] = Quadpoint(+sq15, 0.0, 0.0, 5.0/6.0)
        qps[3] = Quadpoint( +1.0, 0.0, 0.0, 1.0/6.0)

    elif nqp == 5:
        sq37 = math.sqrt(3.0/7.0)

        qps[0] = Quadpoint( -1.0, 0.0, 0.0,  9.0/90.0)
        qps[1] = Quadpoint(-sq37, 0.0, 0.0, 49.0/90.0)
        qps[2] = Quadpoint(  0.0, 0.0, 0.0, 64.0/90.0)
        qps[3] = Quadpoint(+sq37, 0.0, 0.0, 49.0/90.0)
        qps[4] = Quadpoint( +1.0, 0.0, 0.0,  9.0/90.0)


    else:
        print("No lobatto quadrature rule for shp1 and ngp: ", nqp)

    return qps
