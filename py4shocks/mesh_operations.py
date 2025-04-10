# mesh_operations.py
# modified by: Bryan J. Pinargote (bryanpinargote2000@gmail.com)
# mar. 2025

##############################################################################################################

# Import required libraries
import numpy as np

##############################################################################################################

class MeshOperations:


    def __init__(self, mesh):
        self.mesh = mesh


    def create_empty_mesh(self, len_0=1, lenpc_0=1):
        """
        Creates an empty mesh in CGS units and parsec units based on the bounds and dimensions of the input mesh.
        Inputs: len_0 (float, optional): Length unit conversion factor for CGS units. Defaults to 1.
                lenpc_0: Length unit conversion factor for parsec units. Defaults to 1.
        Outputs: tuple: A tuple containing the following arrays:
                   - x, y, z: 1D arrays representing the cell-centered coordinates in CGS units.
                   - xx, yy, zz: 1D arrays representing the vertex coordinates in parsec units.
                   - XX, YY: 2D meshgrid arrays for the x and y coordinates in parsec units.
        """
        # Generate cell-centered coordinates in CGS units
        x  = np.linspace(self.mesh.bounds[0] * len_0, self.mesh.bounds[1] * len_0, self.mesh.dimensions[0] - 1)
        y  = np.linspace(self.mesh.bounds[2] * len_0, self.mesh.bounds[3] * len_0, self.mesh.dimensions[1] - 1)
        z  = np.linspace(self.mesh.bounds[4] * len_0, self.mesh.bounds[5] * len_0, self.mesh.dimensions[2] - 1)
        # Generate vertex coordinates in parsec units
        xx = np.linspace(self.mesh.bounds[0] * lenpc_0, self.mesh.bounds[1] * lenpc_0, (self.mesh.dimensions[0]))
        yy = np.linspace(self.mesh.bounds[2] * lenpc_0, self.mesh.bounds[3] * lenpc_0, (self.mesh.dimensions[1]))
        zz = np.linspace(self.mesh.bounds[4] * lenpc_0, self.mesh.bounds[5] * lenpc_0, (self.mesh.dimensions[2]))
        # Create a 2D meshgrid for x and y coordinates in parsec units
        XX, YY = np.meshgrid(xx, yy)

        return x, y, z, xx, yy, zz, XX, YY
    
    
    def grid_volume(self, lenpc_0=1):
        """
        Computes the volume of each grid cell based on the grid spacing in each direction.
        Inputs: lenpc_0 (float, optional): Length unit conversion factor for parsec units. Defaults to 1.
        Outputs: tuple: A tuple containing:
                    - dxx, dyy, dzz: Grid spacing in the x, y, and z directions, respectively.
                    - dV: Volume of each grid cell.
        """
        # Generate vertex coordinates.
        xx = np.linspace(self.mesh.bounds[0] * lenpc_0, self.mesh.bounds[1] * lenpc_0, (self.mesh.dimensions[0]))
        yy = np.linspace(self.mesh.bounds[2] * lenpc_0, self.mesh.bounds[3] * lenpc_0, (self.mesh.dimensions[1]))
        zz = np.linspace(self.mesh.bounds[4] * lenpc_0, self.mesh.bounds[5] * lenpc_0, (self.mesh.dimensions[2]))

        # Calculate grid spacing in each direction
        dxx = xx[1] - xx[0]       # Grid spacing in the x-direction
        dyy = yy[1] - yy[0]       # Grid spacing in the y-direction
        dzz = zz[1] - zz[0]       # Grid spacing in the z-direction
        
        # Volume
        dV = dxx * dyy * dzz      

        return dxx, dyy, dzz, dV


    def slicing(self, variable_3D, xyz="z", layer=None):
        """
        Extracts a 2D slice from a 3D array at a specified layer along the x, y, or z dimension.
        Inputs: variable_3D (np.ndarray): The 3D numpy array from which the slice will be extracted.
                xyz (str, optional): The dimension along which to slice. Options are "x", "y", or "z". Defaults to "x".
                layer (int, optional): The layer index at which to extract the slice. If None, the middle layer is used. Defaults to None.
        Output: A 2D slice of the input 3D array at the specified layer.
        """
        if layer is None:
            # Calculate the middle layer index if no layer is provided
            if xyz == "x":
                layer = int((self.mesh.dimensions[0] - 1) / 2)  # Middle layer in the x-dimension
            elif xyz == "y":
                layer = int((self.mesh.dimensions[1] - 1) / 2)  # Middle layer in the y-dimension
            elif xyz == "z":
                layer = int((self.mesh.dimensions[2] - 1) / 2)  # Middle layer in the z-dimension

        # Extract the 2D slice based on the specified dimension and layer
        if xyz == "x":
            variable_2D = variable_3D[layer, :, :]  # Slice along the x-dimension
        elif xyz == "y":
            variable_2D = variable_3D[:, layer, :]  # Slice along the y-dimension
        elif xyz == "z":
            variable_2D = variable_3D[:, :, layer]  # Slice along the z-dimension

        return variable_2D