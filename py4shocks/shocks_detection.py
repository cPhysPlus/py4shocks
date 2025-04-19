# shocs_detection.py
# modified by: Bryan J. Pinargote (bryanpinargote2000@gmail.com)
# mar. 2025

##############################################################################################################

# Import required libraries
import os
import yaml
import thermo_fluid_utils
import argparse
import numpy as np
from pyevtk.hl import gridToVTK
from scipy.signal import find_peaks
from mesh_operations import MeshOperations
from read_simulation import SimulationReader

##############################################################################################################

class ShocksCalculations:
    
    def __init__(self, route_vtk_folder, output_vtr_shocks_folder, prefix_name_file, gamma):
        self.output_vtr_shocks_folder = output_vtr_shocks_folder
        self.prefix_name_file = prefix_name_file
        self.route_vtk_folder = route_vtk_folder
        self.gamma = gamma

        # Check if the VTK folder exists and change the working directory to it.
        if os.path.exists(route_vtk_folder):
            os.chdir(route_vtk_folder)
        else:
            raise FileNotFoundError(f"The directory {route_vtk_folder} does not exist.")
        
        # Read normalization units from the simulation log file.
        normalization_units = readunits.read_normalization_units()

        # Extract normalization factors for density, pressure, and velocity.
        self.rho_0 = normalization_units['Density']['number']
        self.prs_0 = normalization_units['Pressure']['number']
        self.v_0 = normalization_units['Velocity']['number']
        len_0 = normalization_units['Length']['number']
        lenpc_0 = len_0 * 3.2408e-19 # Convert cm to parsecs

        # Read the mesh from the initial VTK file.
        mesh_0 = SimulationReader(route_vtk_folder).read_mesh(0, "data")
        empty_mesh = MeshOperations(mesh_0)
        x, y, z, self.xx, self.yy, self.zz, XX, YY = empty_mesh.create_empty_mesh()
        self.dxx, self.dyy, self.dzz, dV = empty_mesh.grid_volume()
        
        # Store the dimensions of the mesh, adjusted for cell-centered data.
        dim = np.asarray(mesh_0.dimensions)
        vec = list(dim)
        
        self.vec = [i-1 for i in dim]


    def min_c(self, css):
        """
        Calculate the minimum sound speed in 3D by comparing neighboring cells.
        Inputs: css (numpy.ndarray): A 3D array of sound speeds.
        Outputs: tuple: A tuple containing the minimum sound speed along each axis (csx, csy, csz).
        """        
        # Move all cells one step to back, left, and down
        cxb1 = np.roll(css, -1, axis=2)  # Roll in the x-direction
        cyl1 = np.roll(css, -1, axis=1)  # Roll in the y-direction
        czd1 = np.roll(css, -1, axis=0)  # Roll in the z-direction
        # Move all cells one step to front, right, and up
        cxf1 = np.roll(css, 1, axis=2)   # Roll in the x-direction
        cyr1 = np.roll(css, 1, axis=1)   # Roll in the y-direction
        czu1 = np.roll(css, 1, axis=0)   # Roll in the z-direction

        # Move all cells two steps to back, left, and down
        cxb2 = np.roll(css, -2, axis=2)  # Roll in the x-direction
        cyl2 = np.roll(css, -2, axis=1)  # Roll in the y-direction
        czd2 = np.roll(css, -2, axis=0)  # Roll in the z-direction
        # Move all cells two steps to front, right, and up
        cxf2 = np.roll(css, 2, axis=2)   # Roll in the x-direction
        cyr2 = np.roll(css, 2, axis=1)   # Roll in the y-direction
        czu2 = np.roll(css, 2, axis=0)   # Roll in the z-direction

        # Find the minimum of the sound speed by comparing rolled arrays
        csxb = np.minimum(cxb1, cxb2)   # Minimum in the x-direction
        csxb = np.minimum(csxb, css)    # Compare with original
        csyl = np.minimum(cyl1, cyl2)   # Minimum in the y-direction
        csyl = np.minimum(csyl, css)    # Compare with original
        cszd = np.minimum(czd1, czd2)   # Minimum in the z-direction
        cszd = np.minimum(cszd, css)    # Compare with original

        csx = np.minimum(cxf1, cxf2)    # Minimum in the x-direction
        csx = np.minimum(csx, csxb)     # Compare with previous minimum
        csy = np.minimum(cyr1, cyr2)    # Minimum in the y-direction
        csy = np.minimum(csy, csyl)     # Compare with previous minimum
        csz = np.minimum(czu1, czu2)    # Minimum in the z-direction
        csz = np.minimum(csz, cszd)     # Compare with previous minimum

        return csx, csy, csz


    def mach_formula(self, v_gr, c, step, gamma):
        """
        Calculate the Mach number using the Rankine-Hugoniot jump conditions that accounts for velocity gradients and specific heat capacity ratio.
        Inputs: v_gr (numpy.ndarray): Velocity gradient array
                c (numpy.ndarray): Local sound speed array
                step (float): Grid spacing (dx, dy, or dz depending on direction)
                gamma (float): Adiabatic index of the gas
        Outputs: m (numpy.ndarray): Array of calculated Mach numbers
        """
        m = abs((-(v_gr * 2. * step * (1. + gamma)) + np.sqrt(16. * c**2 + (v_gr * 2. * step)**2 * (1. + gamma)**2.)) / (4. * c))
        
        return m


    def findpeaks(self, vec, Mx, My, Mz):
        """
        Identify peaks in the directional Mach number arrays.
        Inputs: vec (tuple): A tuple containing the dimensions of the Mach number arrays (depth, height, width).
                Mx (numpy.ndarray): 3D array representing the Mach number in the x-direction.
                My (numpy.ndarray): 3D array representing the Mach number in the y-direction.
                Mz (numpy.ndarray): 3D array representing the Mach number in the z-direction.
        Outputs: A (numpy.ndarray): 3D binary array indicating the presence of peaks in the x-direction.
                 B (numpy.ndarray): 3D binary array indicating the presence of peaks in the y-direction.
                 C (numpy.ndarray): 3D binary array indicating the presence of peaks in the z-direction.
        """

        A = np.zeros((vec[2], vec[1], vec[0]))  # Initialize binary array for x-direction peaks
        B = np.zeros((vec[2], vec[1], vec[0]))  # Initialize binary array for y-direction peaks
        C = np.zeros((vec[2], vec[1], vec[0]))  # Initialize binary array for z-direction peaks
        
        # Find peaks in the y-direction Mach number array
        for i in range(vec[0]):
            for k in range(vec[2]):
                p, _ = find_peaks(My[k, :, i], distance=1)  # Identify peaks
                B[k, p, i] = 1  # Mark peaks in binary array

        # Find peaks in the x-direction Mach number array
        for j in range(vec[1]):
            for k in range(vec[2]):
                p, _ = find_peaks(Mx[k, j, :], distance=1)  # Identify peaks
                A[k, j, p] = 1  # Mark peaks in binary array

        # Find peaks in the z-direction Mach number array
        for i in range(vec[0]):
            for j in range(vec[1]):
                p, _ = find_peaks(Mz[:, j, i], distance=1)  # Identify peaks
                C[p, j, i] = 1  # Mark peaks in binary array

        return A, B, C  # Return binary arrays indicating peak locations


    def save_vtr(self, mach_all, mach_clo, i=0):
        """
        Save Mach number data to a VTR file.
        Inputs: mach_all (numpy.ndarray): 3D array of Mach numbers for all shocks.
                mach_clo (numpy.ndarray): 3D array of Mach numbers for closed shocks.
                i (int, optional): Index for naming the output file. Defaults to 0.
        Outputs: Saves VTR file to disk with the given name (file_name) in the ruote (save_dir).
        """
        # Change to the output directory specified in the class
        os.chdir(self.output_vtr_shocks_folder)
        
        # Create the save directory if it doesn't exist
        if not os.path.exists(self.output_vtr_shocks_folder):
            os.makedirs(self.output_vtr_shocks_folder)
            
        # Construct the full path to save directory
        path = self.output_vtr_shocks_folder
        # Change to the save directory
        os.chdir(path)

        # Save the VTR file:
        # - Construct filename with index i padded to 3 digits
        # - Save mach number arrays as cell data, converting to float32 and transposing
        gridToVTK(self.prefix_name_file + ".0{:03d}".format(i), self.xx, self.yy, self.zz, 
                  cellData={"mach_all": np.float64(mach_all.T), 
                           "mach_clo": np.float64(mach_clo.T)}) # It can be change to float 64 for more pressision


    def shocks_detection(self, i=0):
        """
        Detect shocks in the simulation data for a specific time step.
        Inputs: i (int, optional): Index of the time step. Defaults to 0.
        Outputs: tuple: Arrays of Mach numbers for all shocks and closed shocks.
        """
        sim_object = SimulationReader(self.route_vtk_folder)
        mesh_pluto = sim_object.read_mesh(i, "data")

        # Read simulation variables in CGS units.
        rho_3D_cgs = sim_object.read_variable("rho")
        prs_3D_cgs = sim_object.read_variable("prs")
        vx1_3D_cgs = sim_object.read_variable("vx1")
        vx2_3D_cgs = sim_object.read_variable("vx2")
        vx3_3D_cgs = sim_object.read_variable("vx3")
        tr1_3D = sim_object.read_variable("tr1")
        
        # Calculate sound speed.
        c_so_3D = thermo_fluid_utils.sound_speed(self.gamma, prs_3D_cgs, rho_3D_cgs)
        
        # Find minimum sound speed in each direction.
        csx, csy, csz = self.min_c(c_so_3D)
        
        # Calculate velocity divergence and apply threshold.
        div = thermo_fluid_utils.diver_vel(vx1_3D_cgs, vx2_3D_cgs, vx3_3D_cgs, self.dxx, self.dyy, self.dzz)
        div_tresh = np.where(div[3] >= 0., 0., 1.)
        
        # Calculate pressure gradient modulus and apply threshold.
        mdl_prs_gra = thermo_fluid_utils.pression_gradient_modulus(prs_3D_cgs, self.dzz)
        mdl_prs_gra_tresh = np.where(mdl_prs_gra > 0.1, 1., 0.)
        
        # Apply threshold for tracer variable.
        tr_tresh = np.where(tr1_3D >= 0.10, 1., 0.)
        
        # Combine thresholds to identify shock regions.
        tagtot = div_tresh * mdl_prs_gra_tresh
        tagclo = tagtot * tr_tresh
        
        # Calculate velocity gradients.
        v_gr1 = np.gradient(vx1_3D_cgs, self.dxx, axis=2)
        v_gr2 = np.gradient(vx2_3D_cgs, self.dyy, axis=1)
        v_gr3 = np.gradient(vx3_3D_cgs, self.dzz, axis=0)
        
        # Calculate Mach numbers using the Rankine-Hugoniot jump conditions.
        Mx = self.mach_formula(v_gr1, csx, self.dxx, self.gamma)
        My = self.mach_formula(v_gr2, csy, self.dyy, self.gamma)
        Mz = self.mach_formula(v_gr3, csz, self.dzz, self.gamma)
        
        # Identify peaks in the Mach number arrays.
        Mx1, My1, Mz1 = self.findpeaks(self.vec, Mx, My, Mz)
        
        # Apply peak filtering to Mach numbers.
        Mx *= Mx1
        My *= My1
        Mz *= Mz1
        Ma = np.sqrt(Mx**2 + My**2 + Mz**2) # Calculate the total Mach number.
        
        # Apply thresholds to Mach numbers.
        Mtot = Ma * tagtot
        Mclo = Ma * tagclo
        
        return Mtot, Mclo
    
##############################################################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('--input_route_vtk_folder', type=str, required=False, help="Directory name of the VTK simulation files.")
    parser.add_argument('--output_vtr_shocks_folder', type=str, required=False, help="Directory name to save the VTR file.")
    parser.add_argument('--prefix_name_file', type=str, required=False, help="Base name for the VTR file.")
    parser.add_argument('--gamma', type=float, required=False, help="Isotropic index")

    args = parser.parse_args()
    
    # Use configuration from a YAML file if no command-line arguments are provided
    if args.input_route_vtk_folder == None and args.output_vtr_shocks_folder == None and args.prefix_name_file == None and args.gamma == None:

        with open("config.yaml", "r") as f:
            config = yaml.safe_load(f)

        input_route_vtk_folder = config['route_vtk_folder']
        output_vtr_shocks_folder = config['calculate_shock_cell']['output_vtr_shocks_folder']
        prefix_name_file = config['calculate_shock_cell']['prefix_name_file']
        gamma = config['constants']['gamma']

        readunits = SimulationReader(input_route_vtk_folder)
        shock_detect = ShocksCalculations(input_route_vtk_folder, output_vtr_shocks_folder, prefix_name_file, gamma)
    
    else:
        # Use command-line arguments if provided.
        readunits = SimulationReader(args.input_route_vtk_folder)
        shock_detect = ShocksCalculations(args.input_route_vtk_folder, args.output_vtr_shocks_folder, args.prefix_name_file, args.gamma)
    
    # Process each time step in the simulation.
    for i in range (len(readunits.read_time_steps()[1])):
        print(f"File {i} was saved.")
        Mtot, Mclo = shock_detect.shocks_detection(i)
        shock_detect.save_vtr(Mtot, Mclo, i)
