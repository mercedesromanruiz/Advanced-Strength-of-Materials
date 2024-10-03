# -*- coding: utf-8 -*-
"""
Created on Sat Oct  8 19:59:40 2022

@author: danie
"""

def importGmsh(File):
    D = 3
    file = open(File)
    #D is the dimension of the problem. Vertex's coordinates will return with D values
    mesh = file.read()
    file.close()
    
    
    '''Read the data and split in the different elements to study
    '''
    #Nodes
    N1 = mesh.index('Nodes')
    N2 = mesh.index('EndNodes')
    #In mesh[N1+6] we find the first numeric value just after Nodes
    #In mesh[N2-6] we find the last numeric value just before EndNodes
    Nodes = mesh[N1+6:N2-2]
    Nodes = Nodes +'\n'
    
    #Cells
    C1 = mesh.index('Elements')
    C2 = mesh.index('EndElements')
    #Same idea used with nodes
    Cells = mesh[C1+9:C2-2]
    Cells = Cells +'\n'
    
    #Physical entities
    PN = mesh.count('PhysicalNames')
    if PN == 0:
        print('Physical entities were not found')
    else:
         PN1 = mesh.index('PhysicalNames')
         PN2 = mesh.index('EndPhysicalNames')
         PhysicalNames = mesh[PN1+14:PN2-2]
         PhysicalNames = PhysicalNames + '\n'
         E1 = mesh.index('Entities')
         E2 = mesh.index('EndEntities')
         Entities = mesh[E1+9:E2-2]
         Entities = Entities + '\n'
    
    
    '''Find Physical entities
    '''
    if PN>0:
        #Initilize VerticesEntities and CellsEntities with their dictionary keys
        VerticesEntities = {}
        CellsEntities = {}
        Number_of_lines = PhysicalNames.count('\n')
        End_Line = -1
        #The variable PhysicalNames contains in the first line the total number of PN
        #The next lines contain the dimension of the entity, its tag and its name (each line)
        for i in range(Number_of_lines):
            # For each line we define a vector with the components of that line
            Start_Line = End_Line + 1
            End_Line = PhysicalNames.index('\n',Start_Line)
            Line = PhysicalNames[Start_Line:End_Line]
            L = list(Line.split())
            
            if i == 0:#First line
                DicPN = {}
            else:
                #Include the names of the differents sets in VerticesEntities and CellsEntities
                if int(L[0])==0:
                    VerticesEntities[L[2].replace('"','')] = []
                else:
                    VerticesEntities[L[2].replace('"','')] = []
                    CellsEntities[L[2].replace('"','')] = []
                DicPN.update({int(L[1]):[int(L[0]),L[2]]})
                DicPN[int(L[1])][1] = DicPN[int(L[1])][1].replace('"','')
                #Key: tag of the entity. Contains a list with the dimension of the 
                #entity in the first value, and the name of the PN in the second.
                
        Number_of_lines = Entities.count('\n')
        End_Line = -1
        
        Ver = 0
        Lin = 0
        Sur = 0
        Vol = 0
        #These are auxiliary variables to know how many physical entities are in the
        #Entities section
        Vertex = {}
        Curves = {}
        Surfaces = {}
        Volumes = {}
        #We initialize four dictionaries. Each contains the information of the different
        #entities saved in a different dictionary depending on the dimension
        
        #The variable Entities contains in the first line th number of vertices,
        #lines, surfaces and volumes that are entities. The next lines contain the 
        #information of each entity. Firstly, it includes all the vertices, then the
        #lines, then the surfaces and finally the volumes.
        #Let's talk about the information in each line. For the vertices, the first
        #value of the line is simply to order them, the next three are the coordinates,
        #the next sais if it is a physical entity or not (1 or 0) and then the tag of
        #that PN. For lines, surfaces and volumes, it is the same organization, but
        #instead of three coordinates, it appears 6 geometrical values that are not
        #important in this program. The next values aren't imortant either
        for i in range(Number_of_lines):
            Start_Line = End_Line + 1
            End_Line = Entities.index('\n',Start_Line)
            Line = Entities[Start_Line:End_Line]
            L = list(Line.split())
            
            if i == 0:
                Ver = int(L[0])
                Lin = int(L[1])
                Sur = int(L[2])
                Vol = int(L[3])
            #Firstly, we need to know if the entity is a PN too. We have to evaluate
            #if the forth value of the line is 1
            #Each dictionary (Vertex, Curves, Surfaces and Volumes) contains the information
            #of each entity (the values of the line, in a list), and the key of the 
            #dictionary is the first value of the line
            elif Ver>0:
                if int(L[4])==1:
                    tup = L[1:]
                    Vertex[int(L[0])] = tuple(tup)
                Ver = Ver-1
            elif Lin>0:
                if int(L[7])==1:
                    tup = L[1:]
                    Curves[int(L[0])] = tuple(tup)
                Lin = Lin-1
            elif Sur>0:
                if int(L[7])==1:
                    tup = L[1:]
                    Surfaces[int(L[0])] = tuple(tup)
                Sur = Sur-1
            elif Vol>0:
                if int(L[7])==1:
                    tup = L[1:]
                    Volumes[int(L[0])] = tuple(tup)
                Vol = Vol-1
        
        #At this point we have a dictionary with the different physical entities, and 
        #4 different dictionaries with the vertices, curves, surfaces and volumes that
        #belong to a physical entity
        #This is just saving the data in a more confortable way. After this, when we
        #search the nodes and the cells, we have to include them in the dictionaries
        #VerticesEntities and CellsEntities
    
    '''Find nodes
    '''
    Number_of_lines = Nodes.count('\n')
    End_Line = -1
    Entity = 1
    #This is an auxiliary variable. When an entity is 1-dimensional, the coordinates 
    #of the nodes have 4 components. We have to differenciate whether the line we are
    #studying is an Entity or a Node
    IsPhysical = -1
    #This is an auxiliary variable. When the node belongs to a physical entity, 
    #we have to include them in the VerticesEntities
    Tag = 0
    NNodes = 0
    #This is an auxuliary variable. It counts th number of nodes in the entity when
    #it is 2-dimensional. When Nnodes reaches the number of nodes in the entity, we
    #return the variable Entity the value 1
    tupla = []
    #This is an auxiliary variable. it creates the tuple with the coordinates of 
    #the vertex and the vertex that join a cell, in order to add it to the list
    #of nodes and cells
    Position = []
    #This is an auxiliary variable to know the position of the vertex
    
    #The variable Nodes contains in the first line the number of entities. The
    #number of nodes is the second value. After this, there are three types of
    #lines: 1, three and four-length lines. Th four-length lines contains the
    #information of the entity: dimension, tag, an unuseful value and the number
    #of nodes in the entity. After each four-length line come n lines (n is the 
    #number of nodes in the entity) with the node tag. After all the node tags,
    #the following n lines represent the coordinates of each node.
    #The next loop works as follows. We study each line separately. The first
    #line doesn't have much information. The next is an entity, thats why we initialize
    #the variable entity with 1. If the entity contains nodes, then we have to
    #study those nodes, so entity=0. The lines that contain only the node tag aren't
    #interesting, so we go directly to copy the coordinates. In between, we study
    #if the node belongs to a PN and copy the tag in the variable Tag, while switching
    #the value of IsPhysical to 1. After that, we evaluate if the node studing is
    #a PN and we add it to the VerticesEntities dictionary
    for i in range (Number_of_lines):
        Start_Line=End_Line+1
        End_Line=Nodes.index('\n',Start_Line)
        Line = Nodes[Start_Line:End_Line]
        L = list(Line.split())
        if i==0:
            Nvertices = int(L[1])
            vertices = [None]*Nvertices
             
        elif (Entity==1):
            #The number of nodes in the entity is the last number of the line
            Number_nodes_E=int(L[-1])
            if Number_nodes_E>0:
                Dim = int(L[0])
                if Dim==0:
                    for j in range(len(Vertex)):
                        if Dim==0 and int(L[1])==list(Vertex.keys())[j]:
                            IsPhysical = 0
                            Tag = int(L[1])
                if Dim==1:
                    for j in range(len(Curves)):
                        if int(L[1])==list(Curves.keys())[j]:
                            IsPhysical = 1
                            Tag = int(L[1])
                elif Dim==2:
                    for j in range(len(Surfaces)):
                        if int(L[1])==list(Surfaces.keys())[j]:
                            IsPhysical = 2
                            Tag = int(L[1])
                elif Dim==3:
                    for j in range(len(Volumes)):
                        if int(L[1])==list(Volumes.keys())[j]:
                            IsPhysical = 3
                            Tag = int(L[1])
            
                Entity = 0
                Position = []
        elif (Entity==0):
            if len(L)==1:
                if IsPhysical == 0:
                    VerticesEntities[DicPN[int(Vertex[Tag][4])][1]].append(int(L[0])-1)
                elif IsPhysical == 1:
                    VerticesEntities[DicPN[int(Curves[Tag][7])][1]].append(int(L[0])-1)
                elif IsPhysical == 2:
                    VerticesEntities[DicPN[int(Surfaces[Tag][7])][1]].append(int(L[0])-1)
                elif IsPhysical == 3:
                    VerticesEntities[DicPN[int(Volumes[Tag][7])][1]].append(int(L[0])-1)
                Position.append(int(L[0]))
            else:
                tupla = []
                for j in range(D):  #Write the components of the vertex
                    tupla.append(float(L[j]))
                tupla = tuple(tupla)
                vertices[Position[NNodes]-1] = tupla
                # vertices.append(tupla)
                NNodes = NNodes + 1
                if NNodes>=Number_nodes_E:
                    Entity = 1
                    NNodes = 0
                    IsPhysical = -1
                    
    # #Post processing: prevent, incorrect performance in gmsh
    # x = erase_repeated_nodes(vertices)
    # vertices = x[0]
    #The VerticesEntities dictionary contains lists, but we want tuples
    for i in range(len(VerticesEntities)):
        VerticesEntities[list(VerticesEntities.keys())[i]] = tuple(VerticesEntities[list(VerticesEntities.keys())[i]])
    
    '''Find cells
    '''
    Number_of_lines = Cells.count('\n')
    End_Line = -1
    Entity = 1
    NElem = 0
    cells = []
    
    #The variable Cells is similar to Nodes. The first line contains the number
    #of entities, the total number of elements and two unuseful values.
    #This variable contains also lines containing the information of the entity.
    #This is the same information as before(dimension, tag, unuseful value and 
    #number of elements). After this, the information depends on the dimension of
    #the entity. D=0: number of the vertex, repeated. D=1,2,3: number of the element,
    #number of the two vertices that delimit the element.
    #The next loop works similar to the one for the nodes. We study each line separately.
    #The first line doesn't have much information. The next is an entity, thats
    #why we initialize the variable entity with 1. If the entity contains cells,
    #then we have to study those cells, so entity=0.We need to copy the nodes that
    #delimit the cells. In between, we study if the cell belongs to a PN and copy
    #the tag in the variable Tag, while switching the value of IsPhysical to 1,
    #2 or 3 depending on the dimension of the entity. After that, we evaluate if
    #the cell studing is a PN and we add it to the CellsEntities dictionary
    for i in range (Number_of_lines):
        Start_Line = End_Line+1
        End_Line = Cells.index('\n',Start_Line)
        Line = Cells[Start_Line:End_Line]
        L = list(Line.split())
        
        if i==0: #First line
            Dim = 0
            Ncells = int(L[1])
            cells = [None]*Ncells
        elif (Entity==1):
            #The number of edges in the entity is the last number of the line
            Number_elem_E=int(L[-1])
            if Number_elem_E>0:
                Dim = int(L[0])
                if Dim==1:
                    for j in range(len(Curves)):
                        if int(L[1])==list(Curves.keys())[j]:
                            IsPhysical = 1
                            Tag = int(L[1])
                elif Dim==2:
                    for j in range(len(Surfaces)):
                        if int(L[1])==list(Surfaces.keys())[j]:
                            IsPhysical = 2
                            Tag = int(L[1])
                elif Dim==3:
                    for j in range(len(Volumes)):
                        if int(L[1])==list(Volumes.keys())[j]:
                            IsPhysical = 3
                            Tag = int(L[1])
                Entity = 0
        elif (Entity==0):
            if IsPhysical == 1:
                CellsEntities[DicPN[int(Curves[Tag][7])][1]].append(int(L[0])-1)
            elif IsPhysical == 2:
                CellsEntities[DicPN[int(Surfaces[Tag][7])][1]].append(int(L[0])-1)
            elif IsPhysical == 3:
                CellsEntities[DicPN[int(Volumes[Tag][7])][1]].append(int(L[0])-1)
            if len(L)>1:
                tupla = []
                for j in range(len(L)-1):  #Write the nodes joining the element
                    tupla.append(int(L[j+1])-1)  #Nodes in gmsh start in 1 and in pyris in 0
                tupla = tuple(tupla)
                cells[int(L[0])-1] = tupla
                NElem = NElem + 1
                if NElem>=Number_elem_E:
                    Entity = 1
                    NElem = 0
                    IsPhysical = -1
    
    return(vertices,cells,VerticesEntities,CellsEntities)
