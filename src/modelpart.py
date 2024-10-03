import numpy as np


class class_modelpart():
    """
    This class represents a model part in a finite element analysis. It reads in a mesh file
    and provides information on the mesh geometry, such as the number of cells and points,
    the center and bounds of the mesh, and the physical surface IDs of interest. It also
    sets the IDs and cells of different dimensions and extracts the cells that belong to
    the physical surface.

    Parameters:
    -----------
    mesh_path : str
        The path of the mesh file to be read in.

    Attributes:
    -----------
    mesh : pyvista.PolyData
        The mesh read in from the mesh file.
    cell_dimension : numpy.ndarray
        An array of integers representing the dimension of each cell in the mesh.
    problem_dimension : int
        The dimension of the problem (i.e., the highest dimension of the cells in the mesh).
    vertices : numpy.ndarray
        An array of the vertices of the mesh.
    ids_[n]d : numpy.ndarray (for n = 0,1,2 or 3)
        An array of the IDs of the n-dimensional cells in the mesh.
    cells_[n]d : numpy.ndarray
        An array of the cells with n dimensions in the mesh.
    ids_physical : numpy.ndarray
        An array of the physical surface IDs of interest.

    Methods:
    --------
    get_cell_dimension() -> numpy.ndarray:
        Returns an array of the dimension of each cell in the mesh.
    info() -> None:
        Prints information about the mesh, such as the number of cells and points,
        the number of arrays, the mesh center and bounds.
    set_ids_cells() -> None:
        Sets the IDs and cells of different dimensions based on the dimension of the problem.
    get_physical_surface_cells(number: int) -> numpy.ndarray:
        Returns an array of the cells that belong to the physical surface of interest.
    """

    def __init__(self, mesh_path):

        self.mesh = pv.read(mesh_path)
        self.cell_dimension = self.get_cell_dimension()
        self.problem_dimension = self.cell_dimension[-1]
        self.vertices = self.mesh.points[:,0:self.problem_dimension]

        self.ids_0d = np.array([])
        self.ids_1d = np.array([])
        self.ids_2d = np.array([])
        self.ids_3d = np.array([])
        self.cells_0d =np.array([])
        self.cells_1d =np.array([])
        self.cells_2d =np.array([])
        self.cells_3d =np.array([])

        self.set_ids_cells()


        # Get the physical surface ID of interest
        try:
            self.ids_physical = np.unique(self.mesh['gmsh:physical'])
        except:
            print("No gmsh:physical")

        print(self.problem_dimension)
        #Actual version:
        self.vertices_current = [tuple(sublist) for sublist in self.vertices]

        if self.problem_dimension==1:
            self.cells_current    = [tuple(sublist) for sublist in self.cells_1d]

        elif self.problem_dimension==2:
            self.cells_current    = [tuple(sublist) for sublist in self.cells_2d]

        elif self.problem_dimension==3:
            self.cells_current    = [tuple(sublist) for sublist in self.cells_3d]

    def get_cell_dimension(self):
        cell_dimension = np.zeros(self.mesh.n_cells, dtype=np.int8)
        for i in range(0,self.mesh.n_cells):
            cell_dimension[i] = self.mesh.cell[i].dimension
        return cell_dimension

    def info(self):
        print(f"{'=' * 20} Modelpart {'=' * 20}")
        print(f"Number of Cells: {self.mesh.n_cells}")
        print(f"Number of Points: {self.mesh.n_points}")
        print(f"Number of Arrays: {self.mesh.n_arrays}")
        print(f"Mesh Center: {self.mesh.center}")
        print(f"Mesh Bounds: {self.mesh.bounds}")


    def set_ids_cells(self):

        if self.problem_dimension >= 0:
            self.ids_0d  = np.where(self.cell_dimension == 0)[0]
            self.cells_0d = np.array([self.mesh.cell[i].point_ids for i in self.ids_0d])

        if self.problem_dimension >= 1:
            self.ids_1d  = np.where(self.cell_dimension == 1)[0]
            self.cells_1d = np.array([self.mesh.cell[i].point_ids for i in self.ids_1d])

        if self.problem_dimension >= 2:
            self.ids_2d  = np.where(self.cell_dimension == 2)[0]
            self.cells_2d = np.array([self.mesh.cell[i].point_ids for i in self.ids_2d])

        if self.problem_dimension >= 3:
            self.ids_3d  = np.where(self.cell_dimension == 3)[0]
            self.cells_3d = np.array([self.mesh.cell[i].point_ids for i in self.ids_3d])

    def get_physical_surface_cells(self,number):

        # Find the indices of the cells that belong to the physical surface
        cell_mask = self.mesh['gmsh:physical'] == self.ids_physical[number]
        surf_indices = np.where(cell_mask)[0]

        # Extract the cells that belong to the physical surface
        surf_cells = self.mesh.cells[surf_indices]

        return surf_cells
