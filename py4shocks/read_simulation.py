# read_simulation.py
# modified by: Bryan J. Pinargote (bryanpinargote2000@gmail.com)
# mar. 2025

##############################################################################################################

# Import required libraries
import os
import numpy as np
import pandas as pd
import pyvista as pv

##############################################################################################################

class SimulationReader:


    def __init__(self, route_folder_VTK):
        self.mesh = None
        self.route_folder_VTK = route_folder_VTK


    def read_mesh(self, i, prefix, vtk_or_vtr="vtk"):
        """
        Reads the mesh from a VTK/VTR file located in the specified folder.
        Input: i (int): The index or identifier used to construct the filename. The filename is expected to follow the format "data.0XXX.vtk", where XXX is a zero-padded number based on the input `i`.
        Outputs: mesh: The mesh object read from the VTK/VTR file. The mesh is also stored in the instance variable `self.mesh` for later use.
        """
        # Construct the file path using the provided index `i`
        file_path = os.path.join(self.route_folder_VTK, prefix + ".0{:03d}.".format(i) + vtk_or_vtr)

        # Read the mesh from the VTK/VTR file
        mesh_readed = pv.read(file_path)

        # Store the mesh in the instance variable for future access
        self.mesh = mesh_readed
        
        return self.mesh
    

    def read_variable(self, variable, unit_0=1.):
        """
        Reads the mesh and associated data arrays from a VTK/VTR file and converts them to CGS units.
        Inputs: variable (str): The name of the variable to extract from the mesh data arrays.
                unit_0 (float, optional): Normalization factor to convert values from code units to CGS units.  Defaults to 1.0 (no conversion).
        Outputs: variable_cgs3D (numpy.ndarray): 3D array containing the variable data in CGS units, with dimensions [z_cells, y_cells, x_cells] matching the mesh structure.
        """
        try:
            # Extract the specified variable data array from the mesh and convert to float32 numpy array
            # preference='cell' indicates we want cell-centered data rather than point data
            variable_arr = np.array(pv.get_array(self.mesh, variable, preference='cell'), dtype=np.float32)
            
            # Reshape the 1D array into 3D using the mesh dimensions.
            # Subtract 1 from each dimension because cell data has one less point than vertices.
            variable3D = variable_arr.reshape(self.mesh.dimensions[2] - 1, self.mesh.dimensions[1] - 1, self.mesh.dimensions[0] - 1)
            
            # Convert from code units to CGS units by multiplying by the normalization factor.
            variable_cgs3D = variable3D * unit_0
            
            return variable_cgs3D

        except Exception as e:
            # Handle errors if the mesh is not initialized or other issues occur.
            print('THE MESH IS NOT INITIALIZATE!')
            print('Error: ', e)


    def read_normalization_units(self):
        """
        Reads normalization units from the log file and stores them in a dictionary.
        The log file is expected to contain a section titled "> Normalization Units:" followed by 
        six lines of normalization data. Special handling is applied for the 'Temperature' field to account for additional text.
        Inputs: None
        Outputs: normalization_units (dict): A dictionary where keys are the names of the normalization units (e.g., 'Density', 'Velocity'), and values are dictionaries containing the numerical value and units.
        """
        # Open de .log file
        with open(self.route_folder_VTK + "/Log_Files/pluto.log", 'r') as file:
            normalization_units = {}

            for line in file:
                if line.startswith("> Normalization Units:"):
                    next(file)  # Skip the empty line after the header

                    # Process the next 6 lines containing normalization units
                    for _ in range(6):
                        line = next(file).strip()
                        name, value = line.split(":", 1) # Split only at the first ":"
                        name = name.replace("[", "").replace("]", "").strip()

                        # Special handling for 'Temperature' field.
                        if name == "Temperature":
                            parts = value.strip().split()
                            number = float(parts[0])
                            units = parts[-1].strip("()")
                        else:
                            # Handle other entries.
                            value_part = value.split(",")[0].strip() # Split the primary value (before the comma) into number and unit
                            number_str, units_part = value_part.split(" (", 1)
                            number = float(number_str)
                            units = units_part.replace(")", "").strip()

                        # Add to the dictionary
                        normalization_units[name] = {"number": number, "units": units}

                    break

        return normalization_units
    

    def read_time_steps(self):
        """
        Reads the time column from the file 'vtk.out' and converts the time values to CGS units (seconds) and Mega years (Myr).
        Inputs: None
        Outputs: tuple: A tuple containing two numpy arrays:
                - The first array contains the time values in CGS units (seconds).
                - The second array contains the time values converted to Myr (Myr).
        """
        # Read the 'vtk.out' file
        df_times = pd.read_csv(self.route_folder_VTK + "/vtk.out", sep=r"\s+", header=None)
        
        # Extract the second column, which contains the time values in code units.
        time_code = np.array(df_times.iloc[:, 1])

        # Get the time normalization unit
        tim_0 = self.read_normalization_units()['Time']['number']
        
        # Convert time to CGS units
        tim_cgs = time_code * tim_0
        tim_cgs = np.array(tim_cgs, dtype=[('tim_cgs', 'f8')])  # Add a header 'tim_cgs' to the array.

        # Convert CGS to Myr
        time_myr = tim_cgs['tim_cgs'] / (60. * 60. * 24. * 365.25 * 1e6) # Conversion factors: seconds -> minutes -> hours -> days -> years -> Myr
        time_myr = np.array(time_myr, dtype=[('time_Myr', 'f8')])  # Add a header 'time_Myr' to the array.

        return tim_cgs, time_myr