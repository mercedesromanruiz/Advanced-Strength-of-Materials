#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 13:45:09 2022

@author: ignacio romero
Classes for finite element analysis in pyris
"""

import numpy as np
from numpy import linalg as LA
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import os
from pyevtk.hl import *
from pyevtk.vtk import *


def _openParaviewfile():
    FILE_PATH = "./solution"
    try:
        os.remove(FILE_PATH + ".vtu")
    except:
        pass
    w = VtkFile(FILE_PATH, VtkUnstructuredGrid)
    return w


def _closeParaviewFile(vtkfile):
    vtkfile.save()


def dumpMatrix(filename, a):
    mat = np.matrix(a)
    with open(filename,'wb') as f:
        for line in mat:
            np.savetxt(f, line, fmt='%+.8f')


class Cell:
    """
    Class that encapsules the cell concept, understood as an ordered
    list of vertices.
    """
    def __init__(self, vlist):
        """
        Cell constructor based on an vertex labels. The datum
        vertices contains the list itself; the datum dim gives the
        dimension of the geometry where the cell lives (0-point,
        1-line, 2-surface, 3-volume)
        """
        self.vertices = vlist
        self.dim      = 0


    def __iter__(self):
       """
       Returns the iterator object
       """
       return CellIterator(self)


    def getNVertices(self):
        """
        Number of vertices in a cell
        """
        return len(self.vertices)


    def plot(self, ax):
        """ Plot a cell in the matplotlib format
        """
        p = Polygon(self.vertexCoordinates,
                    closed=True,
                    edgecolor='k',
                    facecolor='coral')
        ax.add_patch(p)


    def print(self):
        """
        Basic information of a cell.
        """
        print("Cell information")
        print("Number of vertices: ", self.getNVertices())
        print("Coodinates of vertices: ")

        #print(self.vertexCoordinates)
        #for k in self.vertexCoordinates:
        #    print(k)


class CellIterator:
    """
    This is an iterator class that allows to loop (with 'for x in cell')
    in a transparent way that proves very convenient.
    """
    def __init__(self, c):
        self._cell = c
        self._index = 0


    def __next__(self):
       """
       Returns the next vertex label.
       """
       if self._index < self._cell.getNVertices() :
           result = self._cell.vertices[self._index]
           self._index += 1
           return result
       # End of Iteration
       raise StopIteration


class Mesh:
    """
    This class embodies the notation of a mesh, understood as a
    geometrical entity that collects all collective information and behavior
    of all its members.
    """

    def __init__(self, vertexCoordinates, cellConnectivities, ns, cs):
        """
        The constructor of the Mesh requires the collection of
        vertices, the cell connectivities, the vertex sets and the cell sets
        """
        self.vertexCoordinates = vertexCoordinates
        self.cells = []
        for c in cellConnectivities:
            self.cells.append(Cell(c))
        self.vertexsets = ns
        self.cellsets = cs


    def cellTypesForParaview(self):
        """
        build the array with cell types that Pareview will employ for plotting
        """
        # Define cell types
        ctype = np.zeros(self.getNCells())
        for c in range(self.getNCells()):
            nv = self.cells[c].getNVertices()
            dim = self.cells[c].dim
            if nv == 1:
                ctype[c] = VtkVertex.tid
            elif nv == 2:
                ctype[c] = VtkLine.tid
            elif nv == 3:
                ctype[c] = VtkTriangle.tid
            elif (nv == 4 and dim == 2):
                ctype[c] = VtkQuad.tid
            elif (nv == 4 and dim == 3):
                ctype[c] = VtkTetra.tid
            elif (nv == 8 and dim == 3):
                ctype[c] = VtkHexahedron.tid
        return ctype


    def dumpMeshForParaview(self, vtkfile):
        """
        This function writes to a file all the information of the mesh
        in vtk format so that post-processing software such as paraview
        --but not only this one-- can use it later.
        """
        x = np.zeros(self.getNVertices())
        y = np.zeros(self.getNVertices())
        z = np.zeros(self.getNVertices())

        for i, c in enumerate(self.vertexCoordinates):
            x[i], y[i] = c

        lv = 0
        conn = []
        offset = np.zeros(self.getNCells())
        for cc, c in enumerate(self.cells):
            conn += c
            lv   += c.getNVertices()
            offset[cc] = lv

        ctype = self.cellTypesForParaview()

        comments = [ "dump mesh for paraview", "pyris 1.0" ]
        IROunstructuredGridToVTK(vtkfile, x, y, z,
                            connectivity = conn,
                            offsets = offset,
                            cell_types = ctype,
                            cellData = None, pointData = None,
                            comments = comments)


    def dumpMeshAndSolutionForParaview(self, solution):
        """ This function dumps in a vtk file all the gometrical information
        of the model as well as the fields obtained from the solution of
        the finite element problem. This file can be used by Paraview --but
        not only this software-- to view the model and its solution.
        """
        FILE_PATH = "./mesh"
        try:
            os.remove(FILE_PATH + ".vtu")
        except:
            pass

        x = np.zeros(self.getNVertices())
        y = np.zeros(self.getNVertices())
        z = np.zeros(self.getNVertices())

        for i, c in enumerate(self.vertexCoordinates):
            x[i], y[i] = c

        lv = 0
        conn = []
        offset = np.zeros(self.getNCells())
        for cc, c in enumerate(self.cells):
            conn += c
            lv   += c.getNVertices()
            offset[cc] = lv

        # Define cell types
        ctype = self.cellTypesForParaview()
        cd = np.random.rand(self.getNCells())
        cellData = {"pressure" : cd}

        if solution[0].size == 1:
            pd = solution
        else:
            vx = np.ascontiguousarray(solution[:,0])
            vy = np.ascontiguousarray(solution[:,1])
            vz = np.zeros(self.getNVertices())
            #pd = np.vstack((vx, vy, vz)).T
            pd = (vx, vy, vz)

        pointData = {"nodal solution" : pd}
        print(pd)

        comments = [ "dump mesh and solution for paraview", "pyris 1.0" ]
        unstructuredGridToVTK(FILE_PATH, x, y, z, connectivity = conn,
                              offsets = offset, cell_types = ctype,
                              cellData = cellData, pointData = pointData,
                              comments = comments)


    def getNVertices(self):
        """ Get the number of vertices in the mesh.
        """
        return len(self.vertexCoordinates)


    def getNCells(self):
        """ Get the number of cells in the mesh.
        """
        return len(self.cells)


    def plot(self):
        """ Plots a simplified mesh """
        self.dumpMeshForParaview()
        #vmesh = pv.read('mesh.vtu')
        #nv = self.getNVertices()
        #vmesh.plot(render_lines_as_tubes=True,
        #           style='wireframe',
        #           line_width=10,
        #           cpos='xy')

        # pl = pv.Plotter()
        # pl.add_mesh(vmesh)

        # labelCoordinates = np.zeros((self.getNVertices(), 3))
        # i = 0
        # for c in self.vertexCoordinates:
        #     labelCoordinates[i,0], labelCoordinates[i,1] = c
        #     i += 1

        # pl.add_point_labels(labelCoordinates, \
        #                        [f'n {i}' for i in range(nv)], \
        #                        font_size=10, point_size=10)

        # pl.camera_position = 'xy'
        # pl.show(render_lines_as_tubes=True)


    def plotDeformed(self, displacement):
        """ Function to plot a deformed mesh. For that the
        arguments must include a (nodal) vector field that
        will be used to translate every vertex.
        """
        self.dumpMeshAndSolutionForParaview(displacement)
        #vmesh = pv.read('mesh.vtu')
        # #nv = self.getNVertices()

        # warped = vmesh.warp_by_vector()
        # pl = pv.Plotter()
        # pl.add_mesh(vmesh)

        # labelCoordinates = np.zeros((self.getNVertices(),3))
        # for i, c in enumerate(self.vertexCoordinates):
        #     labelCoordinates[i,0] = c[0] + displacement[i,0]
        #     labelCoordinates[i,1] = c[1] + displacement[i,1]

        # pl.add_point_labels(labelCoordinates, \
        #                        [f'n {i}' for i in range(nv)], \
        #                        font_size=10, point_size=10)

        # pl.camera_position = 'xy'
        # pl.show()


    def print(self):
        """ Print essential information of the mesh.
        """
        print("\n")
        print("--------------------------------------------------------")
        print("                Mesh information")
        print("--------------------------------------------------------")
        print("Number of vertices    : ", self.getNVertices())
        print("Number of cells       : ", self.getNCells())
        print("Number of vertex sets : ", len(self.vertexsets))
        print("Number of cell sets   : ", len(self.cellsets))


    def printCells(self):
        for c in self.cells:
            c.print()


class Node:
    """The class Node is design to contain a vertex (its position), as
    well as all the finite element solution at this point -- the degrees
    of freedom.
    """

    def __init__(self, label, coordinates):
        """ The constructor of the Node only needs the Vertex that
        it holds.
        """
        self.label  = label
        self.vertex = coordinates
        self.DOFS = []
        self.U = np.empty(shape=(1,0))
        self.Uold = np.empty(shape=(1,0))


    def allocateDOFS(self, el):
        """
        This function allocates the required memory for the
        degrees of freedom of a node that is connected to a
        certain element type.
        """
        dpn = el.theType.getDOFsPerNode()
        #if not self.DOFS:
        #    start = 0
        #else:
        #    start = len(self.DOFS)

        if (dpn > len(self.DOFS)):
            start = len(self.DOFS)
            for n in range(start, dpn):
                self.DOFS.append(-2)

        self.U = np.zeros(len(self.DOFS))
        self.Uold = np.zeros(len(self.DOFS))


    def constrainDOF(self, doflabel, dofvalue):
        """
        Set the degree of freedom with label 'doflabel' to
        be constrained, so it can not changed during the solution.
        This is for boundary conditions.
        """
        self.DOFS[doflabel] = -1
        self.U[doflabel] = dofvalue
        self.Uold[doflabel] = dofvalue


    def coordinates(self):
        """
        Returns the coordinates of the vertex associated with a Node.
        This function should be deleted in the future.
        """
        return self.vertex


    def getNDOFs(self):
        """Get the number of active DOFs in the node structure. It does
        not count constrained DOFs"""
        d = 0
        for n in self.DOFS:
            if (n >= 0):
                d += 1
        return d


    def getNTotalDOFs(self):
        """Get the total number of DOFs in one node structure, counting
        both free and constrained data"""
        return len(self.DOFS)


    def numberConsecutiveDOFs(self, nextdof):
        """After the DOF of a node have been created, and those that
        are constrained have been marked, this function needs to be
        called so that the DOFs are given global numbers, starting
        from the count nextdof.
        """
        for n in range(len(self.DOFS)):
            if (self.DOFS[n] == -2):
                self.DOFS[n] = nextdof
                nextdof += 1
        return nextdof


    def print(self):
        print("Node with label ", self.label, ", coordinates: ", self.vertex)
        print("DOF labels ", self.DOFS, ", value: ", self.U)


class ElementType:
    def __init__(self, familyName):
        self.family = familyName


    def getDOFsPerNode(self):
        print("Careful, parent getDOFsPerNode() should not be called")
        return 2


class Model:

    def __init__(self, aMesh, elmtTypes, cons, loads):
        """
        The Model class holds the Mesh and the data for the
        mechanical analysis. These include the elementTypes,
        the constraints and the loading.
        """
        self.mesh = aMesh
        self.nodes = []
        self.constraints = cons
        self.loading = loads
        self.energy = 0.0
        self.eltypes = elmtTypes
        labeledElmt = {}

        for label, vc in enumerate(aMesh.vertexCoordinates):
            self.nodes.append(Node(label, np.asarray(vc)))

        for setkey in list(aMesh.cellsets.keys()):

            if setkey in elmtTypes:
                cellset = aMesh.cellsets[setkey]
                for cellLabel in cellset:
                    nn = []
                    theCell = aMesh.cells[cellLabel]
                    for n in theCell.vertices:
                        nn.append(self.nodes[n])

                    el = elmtTypes[setkey].createElement(theCell, nn)
                    labeledElmt[cellLabel] = el
                    for n in el.theNodes:
                        n.allocateDOFS(el)

        self.elements = dict(sorted(labeledElmt.items()))
        self.NDOFs = self.numberDOFs()
        self.residual = np.zeros([self.NDOFs])
        self.tangent = np.zeros([self.NDOFs, self.NDOFs])
        self.sparseTangent = lil_matrix((self.NDOFs, self.NDOFs))
        self.Bmatrix = np.zeros([len(self.elements), self.NDOFs])


    def assembleBmatrix(self, B, el, nel):
        """ Assemble the element B-matrix into the global B-matrix.
        This is only valid for truss elements!!!
        The global B matrix has dimension nel x ndofs
        The local B matrix has dimension 1 x 4
        """
        dofs = el.getDOFs()
        for i in range(el.getNDOFs()):
            if (dofs[i] > -1):
                self.Bmatrix[nel ,dofs[i]] += B[0,i]


    def assembleDDEnergy(self, K, el):
        """
        Assemble a matrix obtained from an element into the global
        tangent matrix.
        """
        dofs = el.getDOFs()
        for i in range(el.getNDOFs()):
            if (dofs[i] > -1):
                for j in range(el.getNDOFs()):
                    if (dofs[j] > -1):
                        self.sparseTangent[dofs[i],dofs[j]] += K[i,j]


    def assembleDEnergy(self, elvector, el):
        """
        Assemble a vector obtained from an element into the global
        residual vector.
        """
        dofs = el.getDOFs()
        for d in range(el.getNDOFs()):
            if (dofs[d] > -1):
                self.residual[dofs[d]] += elvector[d]


    def commitCurrentState(self):
        """
        Moves the solution from U to Uold.
        """
        for nd in self.nodes:
            nd.Uold = np.copy(nd.U)


    def dumpSolution(self, step=None):
        try:
            dim = len(self.nodes[0].coordinates())
        except:
            dim = 1

        x = np.zeros(self.mesh.getNVertices())
        y = np.zeros(self.mesh.getNVertices())
        z = np.zeros(self.mesh.getNVertices())

        for i, c in enumerate(self.mesh.vertexCoordinates):
            try:
                dim = len(c)
                x[i]= c[0]
            except:
                dim = 1
                x[i] = c

            #x[i] = c[0]
            y[i] = c[1] if dim>1 else 0.0
            z[i] = c[2] if dim>2 else 0.0

        lv = 0
        conn = []
        offset = np.zeros(self.mesh.getNCells())
        for cc, c in enumerate(self.mesh.cells):
            conn += c
            lv   += c.getNVertices()
            offset[cc] = lv

        # Define cell types
        ctype = self.mesh.cellTypesForParaview()

        # comments = [ "dump mesh for paraview", "pyris 1.0"]
        pointData = {}
        if self.nodes[0].getNTotalDOFs() == 1:
            temp = np.zeros(self.mesh.getNVertices())
            for nc, nd in enumerate(self.nodes):
                temp[nc] = nd.U[0]
            pointData["temperature"] = temp

        else:
            vx = []
            vy = []
            vz = []
            vvel = np.array([np.zeros(3) for i in range(self.mesh.getNVertices())])

            for nc, nd in enumerate(self.nodes):
                vx.append(nd.U[0])
                vy.append(nd.U[1])
                vz.append(0.0)
                vvel[nc] = np.array([nd.U[0], nd.U[1], 0.0])

            pd = ( np.array(vx), np.array(vy), np.array(vz) )
            pointData["displacement"] = pd

        cellData = {}
        resultNames = set()
        for etName in self.eltypes:
            for rname in self.eltypes[etName].resultNames:
                resultNames.add(rname)
        for name in resultNames:
            cd = np.zeros(self.mesh.getNCells())
            for e, elmt in self.elements.items():
                cd[e] = elmt.result(name)
            cellData[name] = cd

        print("Dumping solution to paraview file 'solution.vtu'")
        FILE_PATH = "./solution"
        try:
            os.remove(FILE_PATH + ".vtu")
        except:
            pass

        if step is not None:
            FILE_PATH = "./solution" + str(step)

        unstructuredGridToVTK(FILE_PATH, x, y, z, connectivity = np.array(conn),
                              offsets = offset, cell_types = ctype,
                              cellData = cellData, pointData = pointData)


    def errorNorm(self):
        """ Computes the error in the current solution. In an equilibrium
        point, the energy should be minimum and the error should be zero.
        """
        self.integrateDEnergy()
        return LA.norm(self.residual)


    def getNDOFs(self):
        """ Compute the total number of (free) degrees of freedom in a model.
        """
        d = 0
        for nd in self.nodes:
            d += nd.getNDOFs()
        return d


    def getNAvailableDOFs(self):
        """ the the number of available degrees of freedom in the model. These
        include the true dofs and the constrained dofs """
        d = 0
        for nd in self.mesh.nodes:
            d += nd.getNTotalDOFs()
        return d


    def getNode(self, nodelabel) -> Node:
        """ Return a node in a model, given its label.
        """
        return self.nodes[nodelabel]


    def getSolution(self):
        """ Return an array with the solution in all the nodes """
        c = 0
        sol = np.zeros((self.mesh.getNVertices(), self.nodes[0].getNTotalDOFs()))
        for nd in self.nodes:
            for i in range(nd.U.size):
                sol[c,i] = nd.U[i]
            c += 1
        return sol


    def integrateBmatrix(self):
        """
        Integrates the global B-matrix. Only to be used with trusses, yet
        """
        self.Bmatrix.fill(0.0)
        for k, el in self.elements.items():
            B = el.integrateBmatrix()
            self.assembleBmatrix(B, el, k)


    def integrateEnergy(self) -> float:
        """
        Loops over all elements adding up their internal energy
        and then substract the energy that comes from the external
        loading.
        """
        energy = 0.0
        for k, el in self.elements.items():
            energy += el.integrateEnergy()

        for vertexsetName in self.loading.keys() :
            for forcePair in self.loading[vertexsetName]:
                nodeLabels = self.mesh.vertexsets[vertexsetName]
                for label in nodeLabels:
                    nd = self.getNode(label)
                    u  = nd.U[forcePair[0]]
                    f  = forcePair[1]
                    energy -= f*u

        self.energy = energy
        return energy


    def integrateDEnergy(self):
        """Integrates the gradient of the energy of the whole model, which
        corresponds to the equilibrium equations (i.e. the residual). For
        that it computes the contribution of the elements (internal) and
        then the external loading.
        """
        self.residual.fill(0.0)
        for k, el in self.elements.items():
            DE = el.integrateDEnergy()
            self.assembleDEnergy(DE, el)

        for vertexsetName in self.loading.keys() :
            for forcePair in self.loading[vertexsetName]:
                nodeLabels = self.mesh.vertexsets[vertexsetName]
                for label in nodeLabels:
                    nd  = self.getNode(label)
                    dof = nd.DOFS[forcePair[0]]
                    f   = forcePair[1]
                    if (dof > -1):
                        self.residual[dof] -= f


    def integrateDDEnergy(self):
        """
        Integrates the tangent stiffness matrix by assembling the contributions
        of all the elements in the model.
        """
        self.setTangentToZero()

        for k, el in self.elements.items():
            DDE = el.integrateDDEnergy()
            self.assembleDDEnergy(DDE, el)


    def integrateJet(self, dt) -> float:
        """
        Loops over all elements adding up the transient
        term to the effective energy.
        """
        self.integrateEnergy()

        jet = 0.0
        for k, el in self.elements.items():
            jet += el.integrateJet(dt)

        self.energy += jet
        return self.energy


    def integrateDJet(self, dt):
        """Integrates the gradient of the energy of the whole model, which
        corresponds to the equilibrium equations (i.e. the residual). For
        that it computes the contribution of the elements (internal) and
        then the external loading.
        """
        self.integrateDEnergy()

        for k, el in self.elements.items():
            DE = el.integrateDJet(dt)
            self.assembleDEnergy(DE, el)


    def integrateDDJet(self, dt):
        """Integrates the tangent stiffness matrix by assembling the contributions
        of all the elements in the model.
        """
        self.integrateDDEnergy()

        for k, el in self.elements.items():
            DDE = el.integrateDDJet(dt)
            self.assembleDDEnergy(DDE, el)


    def integrateVolume(self) -> float:
        """
        Loops over all elements adding up their volume
        loading.
        """
        vol = 0.0
        for k, el in self.elements.items():
            vol += el.integrateVolume()

        return vol


    def jetNorm(self, dt):
        """ Computes the error in the current solution. In an equilibrium
        point, the energy should be minimum and the error should be zero.
        """
        self.integrateDJet(dt)
        return LA.norm(self.residual)


    def maximumElementValue(self, name):
        """ Determines the maximum absolute value of the variable 'name' in all
        the elements of the model.
        """
        large = -1e12
        for k, el in self.elements.items():
            large = max(abs(el.result(name)), large)
        return large


    def maximumNodalDOFs(self):
        """
        Determines the maximum absolute value of the degrees of
        freedom in the model, but separated by degree of freedom.
        returns a numpy array of maximum degrees of freedom
        """
        md = np.array([-1.0])
        n = 1
        for nd in self.nodes:
            ndofs = nd.getNTotalDOFs()
            for _ in range(ndofs-n):
                md = np.append(md, -1.0)
                n += 1
            for i in range(ndofs):
                dof = nd.U[i]
                md[i] = max(md[i], abs(dof))

        return md


    def maximumNodalValue(self):
        """
        Determines the maximum absolute value of the degrees of
        freedom in the model.
        """
        return np.amax(np.absolute(self.getSolution()))


    def numberDOFs(self):
        """
        Numbers all the dofs in the model, starting from 0,
        and returning the total number of DOFs
        """
        for vertexsetName in self.constraints.keys() :
            for bcPair in self.constraints[vertexsetName]:
                nodeLabels = self.mesh.vertexsets[vertexsetName]
                for label in nodeLabels:
                    nd = self.getNode(label)
                    nd.constrainDOF(bcPair[0], bcPair[1])

        nextdof = 0
        for k, el in self.elements.items():
            nextdof = el.numberDOFs(nextdof)
        return nextdof


    # def postprocess(self):
    #     pp = pv.read('solution.vtu')
    #     while (var := self.postprocessorQuestion(pp.array_names)) != -1:
    #         pp.plot(scalars = pp.array_names[var],
    #                 component = 0,
    #                 cpos='xy',
    #                 show_scalar_bar = True,
    #                 cmap = 'turbo',
    #                 show_edges = True)


    def postprocessorQuestion(self, vars):
        print("\n\n")
        print("--------------------------------------------------------")
        print("       Select variable for interactive plot")
        print("--------------------------------------------------------\n")

        for i, var in enumerate(vars):
            print("(", i,") --> ",  var)
        print("\nQuit: (q)")

        a = input()
        if a == 'q':
            ia = -1
        else:
            ia = int(a)
        return ia


    def print(self):
        """
        Print essential information of the computational model.
        """
        print("\n")
        print("--------------------------------------------------------")
        print("                Model information")
        print("--------------------------------------------------------")
        print("Number of nodes       : ", len(self.nodes))
        print("Number of elements    : ", len(self.elements))
        print("Number of DOFs        : ", self.getNDOFs())
        # print("Storage of tangent    : ", self.tangent.size * self.tangent.itemsize)
        #print("Storage of sparse     : ", self.sparseTangent.data.nbytes+
        #      self.sparseTangent.indptr.nbytes +
        #      self.sparseTangent.indices.nbytes)

        c = 0
        for vertexsetName in self.constraints.keys() :
            c += len(self.constraints[vertexsetName])
        print("Number of constraints : ", c)
        print("Number of loadings.   : ", len(self.loading))
        print("Total volume          : ", self.integrateVolume())


    def printDetailed(self):
        print("\n")
        print("-----------------------------------------------------------")
        print("               Model detailed information")
        print("-----------------------------------------------------------")
        self.printNodes()
        self.printElements()

        print("\n-----------------------------------------------------------")
        print("                Maximum values")
        print("-----------------------------------------------------------")
        print("Maximum selftress in bars (absolute value): ",
              self.maximumElementValue("sigma"))
        print("Maximum degrees of freedom:", self.maximumNodalDOFs())


    def printElements(self):
        print("Elements:")
        for k, el in self.elements.items():
            print("\nElement ", k)
            el.print()


    def printNodes(self):
        print("Nodes:")
        for nd in self.nodes:
            nd.print()
            print("")


    def printSolution(self):
        c = 0
        for nd in self.nodes:
            print("\nNode ", c, ", DOFS: ")
            c += 1
            for i in nd.U:
                print(i, " ")


    def printEnergy(self):
        """
        Print the value of the potential energy.
        """
        print("Energy      : ", self.energy)


    def printDEnergy(self):
        """
        Print the DEnergy vector.
        """
        print("DEnergy     : ", self.residual)


    def printDDEnergy(self):
        """
        Print the tangent stiffness.
        """
        #print("DDEnergy    : ", self.tangent)
        print("DDEnergy    : ", self.sparseTangent.todense() )


    def setTangentToZero(self):
        """
        Set the entries of the tangent to 0.0
        """
        for k, el in self.elements.items():
            dofs = el.getDOFs()
            for i in range(el.getNDOFs()):
                if (dofs[i] > -1):
                    for j in range(el.getNDOFs()):
                        if (dofs[j] > -1):
                            self.sparseTangent[dofs[i],dofs[j]] = 0.0


class Analysis:

    def __init__(self, aModel):
        self.model = aModel
        self.x = []

    def findLinearizedSolution(self):
        """
        Solves the global system of equations, given the tangent and
        the residual.
        """
        # self.x = np.linalg.solve(self.model.tangent, self.model.residual)
        self.model.sparseTangent = self.model.sparseTangent.tocsr()
        self.x = spsolve(self.model.sparseTangent, self.model.residual)


    def info(self):
        self.model.integrateEnergy()
        self.model.integrateDEnergy()
        print("\n-----------------------------------------------------------")
        print("                 Initial data")
        print("-----------------------------------------------------------")
        print("Energy in the solution: ", self.model.energy)
        print("Error in the solution: ", self.model.errorNorm())


    def localizeSolutionToNodes(self):
        """
        Use the solution self.x to update the degrees of freedom
        in all the dofs of the model.
        """
        for nd in self.model.nodes:
            for i in range(nd.getNTotalDOFs()):
                if (nd.DOFS[i] > -1):
                    nd.U[i] = nd.U[i] - self.x[nd.DOFS[i]]

        self.x.fill(0.0)


    def resetSolution(self):
        for nd in self.model.nodes:
            for i in range(nd.getNTotalDOFs()):
                if (nd.DOFS[i] > -1):
                    nd.U[i] = 0.0

        self.x.fill(0.0)


    def solve(self):
        pass


class StaticLinearAnalysis(Analysis):
    def __init__(self, theModel):
        """ Constructor for the static analysis. It must be called
        with a model, built elsewhere.
        """
        # read the model in the parent class
        super().__init__(theModel)


    def solve(self):
        # solve the linearized problem
        self.model.integrateDEnergy()
        self.model.integrateDDEnergy()
        self.findLinearizedSolution()
        self.localizeSolutionToNodes()

        # Postprocessing
        print("\n-----------------------------------------------------------")
        print("       Solving the linear systems of equations K U = F")
        print("-----------------------------------------------------------")
        self.model.integrateEnergy()
        self.model.integrateDEnergy()
        print("Energy in the solution: ", self.model.energy)
        print("Error in the solution: ", self.model.errorNorm())

        #self.model.mesh.plot()
        #self.model.mesh.plotDeformed(self.model.getSolution())
        #vtkfile = _openParaviewfile()
        #self.model.mesh.dumpMeshForParaview(vtkfile)
        self.model.dumpSolution()
        #_closeParaviewFile(vtkfile)

    def resetSolution(self):
        pass


class TransientAnalysis(Analysis):
    def __init__(self, theModel, dt, tf):
        """ Template for the solution of a linear transient analysis.
        dt : the time step size
        tf : the final time of the integration
        """
        super().__init__(theModel)
        self.dt = dt
        self.tf = tf
        self.model.print()


    def solve(self):
        # information about the intial (zero) solution
        nsteps = int(self.tf//self.dt)
        #self.model.openParaviewfile()
        #self.model.mesh.dumpMeshForParaview()

        for k in range(nsteps):
            print("\n-----------------------------------------------------------")
            print("     Solving step number ", k+1)
            print("-----------------------------------------------------------")
            self.model.integrateJet(self.dt)
            self.model.integrateDJet(self.dt)
            print("Effective energy in the solution: ", self.model.energy)
            print("Error in the transient solution: ", self.model.jetNorm(self.dt))

            # solve the linearized problem
            self.model.integrateDDJet(self.dt)
            self.findLinearizedSolution()
            print("\n...... Solving the linear system of equations K U = F .......\n")
            self.localizeSolutionToNodes()
            self.model.integrateJet(self.dt)
            self.model.integrateDJet(self.dt)
            print("Effective energy in the solution: ", self.model.energy)
            print("Error in the transient solution: ", self.model.jetNorm(self.dt))
            self.model.dumpSolution(k)
            self.model.commitCurrentState()

        # Postprocessing
        #self.model.mesh.plot()
        #self.model.mesh.plotDeformed(self.model.getSolution())
        self.model.printDetailed()
