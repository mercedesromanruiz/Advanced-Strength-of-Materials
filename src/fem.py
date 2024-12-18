#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 13:45:09 2022

@author: ignacio romero
Classes for finite element analysis in pyris
"""

import numpy as np
import math
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
        Cell constructor based on an vertex labels.
        """
        self.vertices = vlist


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


    def print(self):
        """
        Basic information of a cell.
        """
        print("Cell information")
        print("Number of vertices: ", self.getNVertices())
        print("Coodinates of vertices: ")


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


class Constraint:

    def __init__(self,*, dof, value, scaling='t'):
        self.scaling = scaling
        self.dof = dof
        self.value = value


class Loading:

    def __init__(self,*, dof, value, scaling='t'):
        self.scaling = scaling
        self.dof = dof
        self.value = value


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
            dim = len(self.vertexCoordinates[self.cells[c].vertices[0]])
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
        """
        Get the number of vertices in the mesh.
        """
        return len(self.vertexCoordinates)


    def getNCells(self):
        """
        Get the number of cells in the mesh.
        """
        return len(self.cells)


    def plot(self):
        """
        Plots a simplified mesh
        """
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
        """
        Function to plot a deformed mesh. For that the
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
    """
    The class Node is design to contain a vertex (its position), as
    well as all the finite element solution at this point -- the degrees
    of freedom.
    """

    def __init__(self, label, coordinates):
        """
        The constructor of the Node only needs the Vertex that
        it holds.
        """
        self.label  = label
        self.vertex = coordinates
        self.DOFS = []
        self.U = np.empty(shape=(1,0))
        self.Uold = np.empty(shape=(1,0))
        self.V = np.empty(shape=(1,0))
        self.Vold = np.empty(shape=(1,0))


    def allocateDOFS(self, el):
        """
        This function allocates the required memory for the
        degrees of freedom of a node that is connected to a
        certain element type.
        """
        dpn = el.theType.getDOFsPerNode()
        if (dpn > len(self.DOFS)):
            start = len(self.DOFS)
            for n in range(start, dpn):
                self.DOFS.append(-2)

        self.U = np.zeros(len(self.DOFS))
        self.Uold = np.zeros(len(self.DOFS))
        self.V = np.zeros(len(self.DOFS))
        self.Vold = np.zeros(len(self.DOFS))


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
        """
        Get the number of active DOFs in the node structure. It does
        not count constrained DOFs
        """
        d = 0
        for n in self.DOFS:
            if (n >= 0):
                d += 1
        return d


    def getNTotalDOFs(self):
        """
        Get the total number of DOFs in one node structure, counting
        both free and constrained data
        """
        return len(self.DOFS)


    def numberConsecutiveDOFs(self, nextdof):
        """
        After the DOF of a node have been created, and those that
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

    def __init__(self, aMesh, elmtTypes):
        """
        The Model class holds the Mesh and the data for the
        mechanical analysis. These include the elementTypes,
        the constraints and the loading.
        """
        self.mesh = aMesh
        self.nodes = []
        self.constraints = []
        self.loading = []
        self.energy = 0.0
        self.eltypes = elmtTypes
        self.labeledElmt = {}

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
                    self.labeledElmt[cellLabel] = el
                    for n in el.theNodes:
                        n.allocateDOFS(el)


    def addConstraint(self, *, vertexset, dof, value, scaling="t"):
        self.constraints.append((vertexset,
                                 Constraint(dof=dof, value=value, scaling=scaling)))


    def addLoading(self, *, vertexset, dof, value, scaling="t"):
        self.loading.append((vertexset, Loading(dof=dof, value=value, scaling=scaling)))


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
        Moves the solution from tk to tk+1.
        """
        for nd in self.nodes:
            nd.Uold = np.copy(nd.U)
            nd.Vold = np.copy(nd.V)
            nd.V    = np.zeros(nd.getNTotalDOFs())


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

        FILE_PATH = "./solution"
        try:
            os.remove(FILE_PATH + "*.vtu")
        except:
            pass

        if step is not None:
            FILE_PATH = "./solution" + str(step)

        print(f"Dumping solution to paraview file {FILE_PATH}.vtu")
        unstructuredGridToVTK(FILE_PATH, x, y, z, connectivity = np.array(conn),
                              offsets = offset, cell_types = ctype,
                              cellData = cellData, pointData = pointData)


    def errorNorm(self, time):
        """
        Computes the error in the current solution. In an equilibrium
        point, the energy should be minimum and the error should be zero.
        """
        self.integrateDEnergy(time)
        return LA.norm(self.residual)


    def getNDOFs(self):
        """
        Compute the total number of (free) degrees of freedom in a model.
        """
        d = 0
        for nd in self.nodes:
            d += nd.getNDOFs()
        return d


    def getNAvailableDOFs(self):
        """
        Get the number of available degrees of freedom in the model. These
        include the true dofs and the constrained dofs
        """
        d = 0
        for nd in self.mesh.nodes:
            d += nd.getNTotalDOFs()
        return d


    def getNode(self, nodelabel) -> Node:
        """
        Return a node in a model, given its label.
        """
        return self.nodes[nodelabel]


    def getSolution(self):
        """
        Return an array with the solution in all the nodes
        """
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


    def integrateEnergy(self, time) -> float:
        """
        Loops over all elements adding up their internal energy
        and then substract the energy that comes from the external
        loading.
        """
        energy = 0.0
        for k, el in self.elements.items():
            energy += el.integrateEnergy(time)

        for vertexsetName, load in self.loading:
            nodeLabels = self.mesh.vertexsets[vertexsetName]
            t = time
            sc = eval(load.scaling, locals(),
                      {"sin": math.sin, "cos": math.cos, "math": math})
            for label in nodeLabels:
                nd = self.getNode(label)
                u = nd.U[load.dof]
                f = load.value * sc
                energy -= f*u

        self.energy = energy
        return energy


    def integrateDEnergy(self, time):
        """
        Integrates the gradient of the energy of the whole model, which
        corresponds to the equilibrium equations (i.e. the residual). For
        that it computes the contribution of the elements (internal) and
        then the external loading.
        """
        self.residual.fill(0.0)
        for k, el in self.elements.items():
            DE = el.integrateDEnergy(time)
            self.assembleDEnergy(DE, el)

        for vertexsetName, load in self.loading:
            nodeLabels = self.mesh.vertexsets[vertexsetName]
            t = time
            sc = eval(load.scaling, locals(),
                      {"sin": math.sin, "cos": math.cos, "math": math})
            for label in nodeLabels:
                nd = self.getNode(label)
                dof = nd.DOFS[load.dof]
                f = load.value * sc
                if (dof > -1):
                    self.residual[dof] -= f


    def integrateDDEnergy(self, time):
        """
        Integrates the tangent stiffness matrix by assembling the contributions
        of all the elements in the model.
        """
        self.setTangentToZero()

        for k, el in self.elements.items():
            DDE = el.integrateDDEnergy(time)
            self.assembleDDEnergy(DDE, el)


    def integrateJet(self, dt, time) -> float:
        """
        Loops over all elements adding up the transient
        term to the effective energy.
        """
        self.integrateEnergy(time)

        jet = 0.0
        for k, el in self.elements.items():
            jet += el.integrateJet(dt, time)

        self.energy += jet
        return self.energy


    def integrateDJet(self, dt, time):
        """Integrates the gradient of the energy of the whole model, which
        corresponds to the equilibrium equations (i.e. the residual). For
        that it computes the contribution of the elements (internal) and
        then the external loading.
        """
        self.integrateDEnergy(time)

        for k, el in self.elements.items():
            DJ = el.integrateDJet(dt, time)
            self.assembleDEnergy(DJ, el)


    def integrateDDJet(self, dt, time):
        """
        Integrates the tangent stiffness matrix by assembling the contributions
        of all the elements in the model.
        """
        self.integrateDDEnergy(time)

        for k, el in self.elements.items():
            DDJ = el.integrateDDJet(dt, time)
            self.assembleDDEnergy(DDJ, el)


    def integrateVolume(self) -> float:
        """
        Loops over all elements adding up their volume
        loading.
        """
        vol = 0.0
        for k, el in self.elements.items():
            vol += el.integrateVolume()

        return vol


    def jetNorm(self, dt, time):
        """
        Computes the error in the current solution. In an equilibrium
        point, the energy should be minimum and the error should be zero.
        """
        self.integrateDJet(dt, time)
        return LA.norm(self.residual)


    def maximumABSElementValue(self, name):
        """
        Determines the maximum absolute value of the variable 'name' in all
        the elements of the model.
        """
        label = -1
        large = -1e12
        for k, el in self.elements.items():
            if (abs(el.result(name)) > large):
                large = abs(el.result(name))
                label = k
        return large, label


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


    def minimumElementValue(self, name):
        """
        Determines the minimum value of the variable 'name' in all
        the elements of the model and also returns the element label.
        """
        label = -1
        small = float('inf')
        for k, el in self.elements.items():
            et = el.theType
            if (name in et.resultNames):
                if el.result(name) < small:
                    small = el.result(name)
                    label = k

        return small, label


    def numberDOFs(self):
        """
        Numbers all the dofs in the model, starting from 0,
        and returning the total number of DOFs
        """
        for vertexsetName, constraint in self.constraints:
            nodeLabels = self.mesh.vertexsets[vertexsetName]
            for label in nodeLabels:
                nd = self.getNode(label)
                nd.constrainDOF(constraint.dof, constraint.value)

        nextdof = 0
        for k, el in self.elements.items():
            nextdof = el.numberDOFs(nextdof)
        return nextdof


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


    def prepareAnalysis(self):
        self.elements = dict(sorted(self.labeledElmt.items()))
        self.NDOFs = self.numberDOFs()
        self.residual = np.zeros([self.NDOFs])
        self.tangent = np.zeros([self.NDOFs, self.NDOFs])
        self.sparseTangent = lil_matrix((self.NDOFs, self.NDOFs))
        self.Bmatrix = np.zeros([len(self.elements), self.NDOFs])

    def print(self):
        """
        Print essential information of the computational model.
        """
        print("\n")
        print("--------------------------------------------------------")
        print("                Model information")
        print("--------------------------------------------------------")
        print("Number of nodes       : ", len(self.nodes))
        print("Number of eltypes     : ", len(self.eltypes))
        print("Number of DOFs        : ", self.getNDOFs())
        # print("Storage of tangent    : ", self.tangent.size * self.tangent.itemsize)
        #print("Storage of sparse     : ", self.sparseTangent.data.nbytes+
        #      self.sparseTangent.indptr.nbytes +
        #      self.sparseTangent.indices.nbytes)
        print("Number of constraints : ", len(self.constraints))
        print("Number of loadings    : ", len(self.loading))
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

        mm, label = self.maximumABSElementValue("sigma")
        print(f"Maximum stress (absolute value): {mm} (element {label})")
        print("Maximum degrees of freedom:", self.maximumNodalDOFs())
        safety, label = self.minimumElementValue("safety")
        print(f"Minimum safety factor: {safety} (element {label})")


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
        self.model.prepareAnalysis()
        self.x = []


    def findLinearizedSolution(self):
        """
        Solves the global system of equations, given the tangent and
        the residual.
        """
        # self.x = np.linalg.solve(self.model.tangent, self.model.residual)
        self.model.sparseTangent = self.model.sparseTangent.tocsr()
        self.x = spsolve(self.model.sparseTangent, self.model.residual)


    def info(self, time):
        self.model.integrateEnergy(time)
        self.model.integrateDEnergy(time)
        print("\n-----------------------------------------------------------")
        print("    Info at time t:", time)
        print("-----------------------------------------------------------")
        print("Energy in the solution: ", self.model.energy)
        print("Error in the solution: ", self.model.errorNorm())


    def localizeSolutionToNodes(self, dt):
        """
        Use the solution self.x to update the degrees of freedom
        in all the dofs of the model.
        """
        for nd in self.model.nodes:
            for i in range(nd.getNTotalDOFs()):
                globalEq = nd.DOFS[i]
                if (globalEq > -1):
                    nd.U[i] = nd.U[i] - self.x[globalEq]
                    nd.V[i] = (nd.U[i] - nd.Uold[i])/dt

        self.x.fill(0.0)


    def resetSolution(self):
        for nd in self.model.nodes:
            for i in range(nd.getNTotalDOFs()):
                if (nd.DOFS[i] > -1):
                    nd.U[i] = 0.0
                    nd.V[i] = 0.0

        self.x.fill(0.0)


    def solve(self):
        pass


class StaticLinearAnalysis(Analysis):
    def __init__(self, theModel, dt=1.0, tf=1.0):
        """
        Constructor for the static analysis. It must be called
        with a model, built elsewhere.
        """
        super().__init__(theModel)
        self.dt = dt
        self.tf = tf


    def solve(self, label=None):
        nsteps = int(self.tf//self.dt)
        t = 0.0
        self.model.dumpSolution(0)
        for k in range(nsteps):
            t += self.dt

            if (label):
                ss = label
            else:
                ss = k+1

            print("\n-----------------------------------------------------------")
            print("     Static analysis. Solving step number", ss, ", time:", t)
            print("-----------------------------------------------------------")
            self.model.integrateEnergy(t)
            self.model.integrateDEnergy(t)
            self.model.integrateDDEnergy(t)

            print("Effective energy in the solution: ", self.model.energy)
            print("Error in the solution: ", self.model.errorNorm(t))

            self.findLinearizedSolution()
            self.localizeSolutionToNodes(t)
            self.model.integrateEnergy(t)
            self.model.integrateDEnergy(t)
            print("Energy in the solution: ", self.model.energy)
            print("Error in the solution: ", self.model.errorNorm(t))
            self.model.dumpSolution(step=ss)


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
        nsteps = int(self.tf//self.dt)
        t = 0.0
        self.model.dumpSolution(0)
        for k in range(nsteps):
            t += self.dt
            print("\n-----------------------------------------------------------")
            print("     Solving step number", k+1, ", time:", t)
            print("-----------------------------------------------------------")
            self.model.integrateJet(self.dt, t)
            self.model.integrateDJet(self.dt, t)
            print("Effective energy in the solution: ", self.model.energy)
            print("Error in the transient solution: ",
                  self.model.jetNorm(self.dt, t))

            # solve the linearized problem
            self.model.integrateDDJet(self.dt, t)
            self.findLinearizedSolution()
            print("\n...... Solving the linear system of equations K U = F .......\n")
            self.localizeSolutionToNodes(self.dt)
            self.model.integrateJet(self.dt, t)
            self.model.integrateDJet(self.dt, t)
            print("Effective energy in the solution: ", self.model.energy)
            print("Error in the transient solution: ",
                  self.model.jetNorm(self.dt, t))
            self.model.dumpSolution(k+1)
            self.model.commitCurrentState()


class EnergySteppingAnalysis(Analysis):
    def __init__(self, theModel, dt, tf):
        """ Template for the solution of a linear transient analysis.
        dt : the time step size
        tf : the final time of the integration
        """
        super().__init__(theModel)
        self.dt = dt
        self.tf = tf
        self.model.print()
        self.TOL = 1e-5
        self.OMT = 1.0-TOL

    def energySteppingFindT(q0, p0, L0, h, target,
                            rtol=1e-6, xtol=1e-5, maxiter=15, debug=False):
        """
        Find the time at which the system crosses the next level set of the
        terraced potential energy when the system moves in the direction of the
        velocity.

        The function uses the secant method when it detects crossing a terraced level
        in the direction of the gradient. In a basin, it switches to brendt method.

        q0, p0: initial position and momentum
        L0: initial level of terraced potential (need not be equal to V(q0) !!)
        h:  terraced potential steps
        target: target distribution. A function that must have two methods implemented
                target.potential(q) and target.gradPotential(q)
        rtol: relative tolerace in the function evaluation for convergence of root
        xtol: relative tolerance in the unknown for convergence of root finding
        maxiter: maximum number of iterations in root finding

        """
        V0 = L0*h
        W0 = target.potential(q0)
        g0 = target.gradPotential(q0)
        gradNorm = np.linalg.norm(g0)
        v0 = p0
        v0norm = np.linalg.norm(v0)
        a0 = np.dot(v0, g0)

        if debug:
            print("    Find T, L0: {:d}, L(q0): {:.4e}, angle: {:.4e}".format(
                L0, W0/h, a0))

        # find bracket [t0, t] where energy step is met. Determine up/low bounds
        Gtop = lambda z: target.potential(q0+z*v0) - V0 - h
        Gbot = lambda z: target.potential(q0+z*v0) - V0

        dx = min(h/gradNorm,1.0)
        t = 1.0e-3*dx/v0norm

        # choose t0 to be slightly positive, but make sure that
        # the starting point is good ... Gtop<0 and Gbot>0
        t0 = t
        giter = 0
        maxgiter = 10
        while ((Gtop(t0)>0 or Gbot(t0)<0) and giter<maxgiter):
            t0 *= 0.5
            giter += 1

        if giter == maxgiter:
            print("\n Number of splits for initial t is too large.", giter)
            print(" Gtop(0): {:.4e}, Gbot(0): {:.4e}".format(Gtop(0), Gbot(0)))
            sys.exit(0)

        # determine whether the trajectory crosses up or downward energy
        direction = 0
        maxgiter = 10*maxiter
        giter = 0

        Gtop0 = Gtop(t0)
        Gbot0 = Gbot(t0)

        while (Gtop(t)*Gtop0>0 and Gbot(t)*Gbot0>0 and giter<maxgiter):
            t *= 1.61
            giter += 1

        if giter == maxgiter:
            print("\n Number of iterations in up/down detection is too large.", giter)
            sys.exit(0)
        elif (Gtop(t)*Gtop0<0):
            direction = +1
            G = Gtop
        else:
            direction = -1
            G = Gbot

        if debug and direction==1:
            print("    Gtop[{:.4e},{:,.4e}] = [{:.4e},{:,.4e}], dir={:,.1g}".format(
                        0, t, Gtop(0), Gtop(t), direction))
            print("    Gtop[{:.4e},{:,.4e}] = [{:.4e},{:,.4e}], dir={:,.1g}".format(
                        t0, t, Gtop(t0), Gtop(t), direction))

        if debug and direction==-1:
            print("    Gbot[{:.4e},{:,.4e}] = [{:.4e},{:,.4e}], dir={:,.1g}".format(
                        0, t, Gbot(0), Gbot(t), direction))
            print("    Gbot[{:.4e},{:,.4e}] = [{:.4e},{:,.4e}], dir={:,.1g}".format(
                        t0, t, Gbot(t0), Gbot(t), direction))

        # backup initial point to ensure sign change is in bracket
        #t0 = -direction*mysign(a0)*1e-3*dx/v0norm

        if (G(t0)*G(t)>=0):
            print("\n*Gtop[{:.4e},{:,.4e}] = [{:.4e},{:,.4e}], dir={:,.1g}".format(
                            t0, t, Gtop(t0), Gtop(t), direction))
            print("*Gbot[{:.4e},{:,.4e}] = [{:.4e},{:,.4e}], dir={:,.1g}".format(
                            t0, t, Gbot(t0), Gbot(t), direction))
            print("*G   [{:.4e},{:,.4e}] = [{:.4e},{:,.4e}], dir={:,.1g}".format(
                            t0, t, G(t0), G(t), direction))

        secant = False
        if secant:
            sol = root_scalar(G, x0=t0, x1=t, method="secant",
                              rtol=rtol, xtol=xtol*dx, maxiter=maxiter)
        else:
            sol = root_scalar(G, bracket=[t0,t], method="brentq",
                              rtol=rtol, maxiter=maxiter)

        if (not(sol.converged) and debug):
            print("\nFirst trial error in root finding for energy-stepping update")
            print("Cause of failure:", sol.flag)

        if debug:
            print(" -> Time increment found {0:.4f} ({1:d} iterations), {2:s}".format(
                sol.root, sol.iterations, sol.flag))

        if (sol.root < 1e-7):
            print("    dt {:.4e} is negative or suspiciously small...".format(sol.root))
            print("    Find T, L0: {:.4e}, L(q0): {:.4e}, Grad: {:.4e}".format(
                L0, W0/h, gradNorm))
            sys.exit(0)

        return sol.root


    def energySteppingProposal(q0, p0, h, T, target, debug=False):
        """
        Perform integration for a period of time of length T of a Hamiltonian system
        usign the energy-stepping method. In this method, the total energy
        is terraced in steps of height h.

        input:
          - q0, p0 : initial phase coordinates
          - h: terrace height
          - T: integrating time in each step
          - target: target probability distribution

        returns:
          - position q, momentum p
          - terraced potential level at time T
          - sliding motions
          - the total time integrated

        notes:
           - the potential energy need not be the potential energy at q due
             to the terraced formulation
           - sliding motions: movements that changed q (not simply
             climbed or dropped a step. They require a gradV evaluation)

           Vn: terraced potential at qn = h * floor(V(qn))
           Ln: steps of the terraced potential at qn ... floor(Vn/h)
           Wn: potential at qn = V(qn)
           angle an: angle between the velocity vn and the gradient gn
        """
        tf = 0.0
        L0 = math.floor(target.potential(q0)/h)
        V0 = L0*h
        K0 = 0.5*np.dot(p0, p0)

        qn = np.copy(q0)
        vn = np.copy(p0)
        Ln = L0
        Vn = V0
        case = 0
        slided = 0
        step = 1

        while tf<T:

            # note the difference between Vn, the terraced energy and
            # Wn, the true energy at qn
            gn = target.gradPotential(qn)
            gradnorm = np.linalg.norm(gn)
            nn = gn/gradnorm
            an = np.dot(vn, nn)
            Wn = target.potential(qn)
            dt = 0.0

            if debug:
                print("\n>>> Integrate time step {:d}".format(step))
                print("    Initial data. (h: {:.4e})".format(h))
                print("    Vn: {:.4e}, Ln: {:d}, Wn: {:.4e}, L(qn): {:.4e}, angle: {:.4e}".format(
                    Vn, Ln, Wn, Wn/h, an))
                print("    Test1, (Wn-Vn)/h: {:.4e}  - {:}".format((Wn-Vn)/h, Wn-Vn>= OMT*h))
                print("    Test2, an: {:.4e}         - {:}".format(an, an>0))
                print("    Test3, an*an/(2*h): {:.4e}- {:}".format(an*an/(2*h), an*an>2*h))

            # compute an intermediate position qh and velocity vh, placing the point on a terrace
            # with the velocity away from any wall
            # point is facing wall and has enough energy to go up
            if (Wn-Vn >= OMT*h and an>0.0 and an*an >= 2.0*h):
                vh = vn + (-an + math.sqrt(an*an-2.0*h))*nn
                Lh = Ln + 1
                if debug: print("    Climbed wall. Initial level: {:d}".format(Lh))

            # point is facing wall, but has no energy to go up
            elif (Wn-Vn >= OMT*h and an>0.0 and an*an<2.0*h):
                vh = vn - 2.0*an*nn
                Lh = Ln
                if debug: print("    Bumped wall. Initial level: {:d}".format(Lh))

            # point is facing cliff, so go down
            elif (Wn-Vn < TOL*h and an<0.0):
                vh = vn + (-an - math.sqrt(an*an+2.0*h))*nn
                Lh = Ln - 1
                if debug: print("    Fell cliff. Initial level: {:d}".format(Lh))

            # point is pointing in the right direction
            else:
                vh = np.copy(vn)
                Lh = Ln
                if debug: print("    Moving ok. Initial level: {:d}".format(Lh))

            # point is in the flat area and with the right direction
            # Move until hit wall or fall cliff
            slided += 1
            qh = np.copy(qn)
            rtol = h*TOL/max(1000,gradnorm)
            dt = energySteppingFindT(qh, vh, Lh, h, target, rtol=rtol, debug=debug)

            if (tf+dt > T):
                dt = T-tf
                if debug: print("    Slide and stop (dt={:.4e})".format(dt), end="")

            q1 = qh + vh*dt
            v1 = np.copy(vh)
            L1 = Lh
            V1 = Lh*h

            Kn = 0.5*np.dot(vn,vn)
            En = Vn + Kn

            K1 = 0.5*np.dot(v1,v1)
            E1 = V1 + K1

            if debug:
                print("    End of update with dt: {:5g}".format(dt))
                print("    Ln: {:d}, L(qn): {:.4g}, Kn: {:.4e}".format(
                    Ln, target.potential(qn)/h, Kn))
                print("    L1: {:d}, L(q1): {:.4g}, K1: {:.4e}".format(
                    L1, target.potential(q1)/h, K1))

            if (math.fabs(En-E1)>1e-4*h):
                print("Error in the energy update! ", E1-En, ", case:", case)
                print("--")
                sys.exit(5)

            qn = np.copy(q1)
            vn = np.copy(v1)
            Vn = V1
            Ln = L1
            tf += dt
            step += 1

        if debug:
            print("\n*** Proposal is computed")
            print("    L1: {:.4e}, L(q1): {:.4e}, K1: {:.4e}".format(
                L1, target.potential(q1)/h, K1))
            print("    H0: {:.4e}, H1: {:.4e}".format(V0+K0, L1*h+K1))

        return q1, v1, L1, slided, tf


    def energySteppingMC(target, initial, iterations=5000,
                         h=0.1, T=1.0, debug=False):
        """
        Generate a random sample distributed according to a
        given probability distribution using the Energy-Stepping MC method.

        Inputs:
          - target: the distribution to sample from.
          - initial: generator of the first sample.
          - iterations: number of samples in the Markov chain.
          - h: potential energy terrace size
          - T: time to integrate in each proposal
          - debug: optional debug for print information

        Returns:
          - the chain
          - the ratio of slided motions / MC iterations
          - the number of gradV evaluations
        """
        samples = initial
        slided_total = 0
        dim = initial.shape[1]
        eye = np.identity(dim)
        zz = np.zeros(dim)

        for it in progressbar(range(iterations)):
            qn = samples[-1]
            pn = np.random.multivariate_normal(mean=zz, cov=eye)

            if debug:
                print("\n===========================================")
                print("              New proposal ({:d})".format(it))
                print("===========================================")

            qn1, pn1, Ln1, slided, tf = energySteppingProposal(qn, pn, h, T,
                                                               target, debug)

            samples = np.vstack((samples, qn1))
            slided_total += slided

        return samples, slided_total/iterations, slided_total


        def solve(self):
            nsteps = int(self.tf//self.dt)
            t = 0.0
            self.model.dumpSolution(0)
            for k in range(nsteps):
                t += self.dt
                print("\n-----------------------------------------------------------")
                print("     Solving step number", k+1, ", time:", t)
                print("-----------------------------------------------------------")
                self.model.integrateJet(self.dt, t)
                self.model.integrateDJet(self.dt, t)
                print("Effective energy in the solution: ", self.model.energy)
                print("Error in the transient solution: ",
                      self.model.jetNorm(self.dt, t))

                # solve the linearized problem
                self.model.integrateDDJet(self.dt, t)
                self.findLinearizedSolution()
                print("\n...... Solving the linear system of equations K U = F .......\n")
                self.localizeSolutionToNodes(self.dt)
                self.model.integrateJet(self.dt, t)
                self.model.integrateDJet(self.dt, t)
                print("Effective energy in the solution: ", self.model.energy)
                print("Error in the transient solution: ",
                      self.model.jetNorm(self.dt, t))
                self.model.dumpSolution(k+1)
                self.model.commitCurrentState()
