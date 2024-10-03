#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 21:04:09 2022

@author: ignacio romero
"""

def processEntities(entitiesLines, debug = False):
    """ obtain the points, curves, surfaces and volumes
    """
    points = {}
    l = entitiesLines[0].split()
    np, nc, ns, nv = int(l[0]), int(l[1]), int(l[2]), int(l[3])

    # read points in a dictionary points[tag] = physical entities of this point
    k = 1
    for a in range(np):
        l = entitiesLines[k].split()
        tag = int(l[0])
        #x, y, z = float(l[1]), float(l[2]), float(l[3])
        nphys = int(l[4])
        physTags = []
        for b in range(nphys):
            physTags.append( int(l[5+b]) )
        points[tag] = physTags
        k += 1

    # read curves in a dictionary points[tag] = physical entities of this line
    curves = {}
    for a in range(nc):
        l = entitiesLines[k].split()
        tag = int(l[0])
        #x, y, z = float(l[1]), float(l[2]), float(l[3])
        #X, Y, Z = float(l[4]), float(l[5]), float(l[6])
        nphys = int(l[7])
        physTags = []
        for b in range(nphys):
            physTags.append( int(l[8+b]) )
        curves[tag] = physTags
        k += 1

    # read surfaces in a dictionary surfaces[tag] = physical entities of this surface
    surfaces = {}
    for a in range(ns):
        l = entitiesLines[k].split()
        tag = int(l[0])
        #x, y, z = float(l[1]), float(l[2]), float(l[3])
        #X, Y, Z = float(l[4]), float(l[5]), float(l[6])
        nphys = int(l[7])
        physTags = []
        for b in range(nphys):
            physTags.append( int(l[8+b]) )
        surfaces[tag] = physTags
        k += 1

    # read volumes in a dictionary volumes[tag] = physical entities of this surface
    volumes = {}
    for a in range(nv):
        l = entitiesLines[k].split()
        tag = int(l[0])
        #x, y, z = float(l[1]), float(l[2]), float(l[3])
        #X, Y, Z = float(l[4]), float(l[5]), float(l[6])
        nphys = int(l[7])
        physTags = []
        for b in range(nphys):
            physTags.append( int(l[8+b]) )
        volumes[tag] = physTags
        k += 1

    if debug:
        print("Points: ", points)
        print("Curves: ", curves)
        print("Surfaces: ", surfaces)
        print("Volumes: ", volumes)
        print()

    return points, curves, surfaces, volumes



def processNodes(nodesLines, surfaces, volumes, vsets, tag2PN, gDim, debug = False):
    """Based on the information on the nodesLines, and extra ,
    build the vertices dictionary
    """
    vertices = {}
    l = nodesLines[0].split()
    nb, nn = int(l[0]), int(l[1])
    #mt, Mt = int(l[2]), int(l[3])

    k = 1
    for a in range(nb):
        l = nodesLines[k].split()
        k += 1
        entityDim, entityTag = int(l[0]), int(l[1])
        #parametric = int(l[2])
        nn = int(l[3])
        ntags = []
        for b in range(nn):
            l = nodesLines[k].split()
            k += 1
            nodetag = int(l[0])
            ntags.append(nodetag)

            if entityDim == 2 :
                physTags = surfaces[entityTag]
                for e in physTags:
                    vsets[ tag2PN[e] ] += (nodetag,)
            if entityDim == 3 :
                physTags = volumes[entityTag]
                for e in physTags:
                    vsets[ tag2PN[e] ] += (nodetag,)
        for b in range(nn):
            l = nodesLines[k].split()
            k += 1
            if gDim == 1: vertices[ntags[b]] = (float(l[0]))
            if gDim == 2: vertices[ntags[b]] = (float(l[0]), float(l[1]))
            if gDim == 3: vertices[ntags[b]] = (float(l[0]), float(l[1]), float(l[2]))

    if debug:
         print("Vsets: ", vsets)
         print()

    return vertices



def processElements(elementsLines, surfaces, volumes, csets, tag2PN, gDim, debug=False):
    cells = {}
    l = elementsLines[0].split()
    nb, ne = int(l[0]), int(l[1])
    #mt, Mt = int(l[2]), int(l[3])

    k = 1
    for a in range(nb):
        l = elementsLines[k].split()
        k += 1
        entityDim, entityTag, eltype, ne = int(l[0]), int(l[1]), int(l[2]), int(l[3])
        nodetags = []
        for b in range(ne):
            l = elementsLines[k].split()
            k += 1
            tag = int(l[0])
            nodetags = (int(l[1]),)
            if   eltype == 1: nn = 2
            elif eltype == 2: nn = 3
            elif eltype == 3: nn = 4
            elif eltype == 4: nn = 4
            elif eltype == 15: nn = 1
            else:
                print("element type not defined for ", eltype)
            for c in range(nn-1):
                nodetags = nodetags + (int(l[2+c]),)
            cells[tag] = nodetags

            if entityDim == 2 :
                physTags = surfaces[entityTag]
                for e in physTags:
                    csets[ tag2PN[e] ] += (tag,)
            if entityDim == 3 :
                physTags = volumes[entityTag]
                for e in physTags:
                    csets[ tag2PN[e] ] += (tag,)

    if debug:
        print("Cells: ", cells)
        print()
        print("Cell sets", csets)
        print()

    return cells




def remap(vertices, cells, vsets, csets, debug = False):
    """Remap gmsh data so that there are no gaps in the numbering
    """

    #count the number of vertices
    unique_vertices = set()
    for v in vertices:
        unique_vertices.update(vertices.keys())

    rvertices = []
    rvsets = []
    if len(unique_vertices) == len(vertices):
        rvertices = vertices
        rvsets = csets

    unique_cells = set()
    for c in cells:
        unique_cells.update(cells.keys())

    rcells = []
    rcsets = []
    if len(unique_cells) == len(cells):
        rcells = cells
        rcsets = csets

    if debug and False:
        print(unique_vertices)
        print("Number of unique v", len(unique_vertices))

        print(unique_cells)
        print("Number of unique c", len(unique_cells))

    return rvertices, rcells, rvsets, rcsets




def gmshImport(filename, debug=False, dim=3):
    """Read a file in gmsh format, building the arrays for the vertices, the
    cells, the vertex sets, and the cell sets.
    Arguments are:
        thefile, typically with a .msh extension
        debug, set equal to True to display extra messages
        dim, dimension of the body (default =3).
    """
    gDim = dim

    # read all the data as lists of ints and strings
    physNames = fileBlock(filename, "$PhysicalNames", "$EndPhysicalNames")
    entitiesLines = fileBlock(filename, "$Entities", "$EndEntities")
    nodesLines    = fileBlock(filename, "$Nodes", "$EndNodes")
    elementsLines = fileBlock(filename, "$Elements", "$EndElements")

    if debug and False:
        print("Physical names: ", physNames)
        print("Entities lines: ", entitiesLines)
        print("Nodes lines: ", nodesLines)
        print("Elements lines: ", elementsLines)

    # extract physical names as two maps, name->dimension, tag->name
    # initialize vertex sets and cell sets
    PNdim = {}
    tag2PN = {}
    vsets = {}
    csets = {}
    for f in physNames[1:]:
        l = f.split()
        dim, tag = int(l[0]),  int(l[1])
        name = l[2][1:-1]
        PNdim[name] = dim
        tag2PN[tag] = name
        vsets[name] = ()
        csets[name] = ()

    if debug and False:
        print("Physical names dimensions: ", PNdim)
        print("Physical namnes tags: ", tag2PN)


    points, curves, surfaces, volumes = processEntities(entitiesLines, debug)
    vertices = processNodes(nodesLines, surfaces, volumes, vsets, tag2PN, gDim, debug)
    cells = processElements(elementsLines, surfaces, volumes, csets, tag2PN, gDim, debug)

    # change numbering, if necessary, to avoid gaps
    newvertices, newcells, newvsets, newcsets = remap(vertices, cells, vsets, csets, debug)

    return newvertices, newcells, newvsets, newcsets



def fileBlock(filename, startstr, endstr):
    """ return a list of lines in the file from the line starting with the
    string startstr, up to the line ending with the string endstr. It does
    not return the first and the last """
    file = open(filename, 'r')

    block = []
    for line in file:
        if line != None and startstr in line:
            line = file.readline()
            while line != None and not(endstr in line):
                block.append(line.strip())
                line = file.readline()

    file.close()
    return block


gmshImport("examples/test.msh", debug = False)
