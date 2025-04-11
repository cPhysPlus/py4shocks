# shocs_detection.py
# modified by: Bryan J. Pinargote (bryanpinargote2000@gmail.com)
# mar. 2025

##############################################################################################################

# Import required libraries
import yaml
import argparse
import numpy as np
import pandas as pd
from mesh_operations import MeshOperations
from read_simulation import SimulationReader

##############################################################################################################

class last_index_position:


    def __init__(self, input_route_vtk_folder, input_folder_mach, name_tracking_output, output_tracking_dir, name_mach_file, threshold_tracer, threshold_rho, threshold_mach):
        self.input_route_vtk_folder = input_route_vtk_folder
        self.input_folder_mach = input_folder_mach
        self.name_mach_file = name_mach_file
        self.output_tracking_dir = output_tracking_dir
        self.name_tracking_output = name_tracking_output
        self.threshold_tracer = threshold_tracer
        self.threshold_rho = threshold_rho
        self.threshold_mach = threshold_mach
    

    def distance_travel_2D(self, arr_2D, mesh):
        """
        Function to calculate the maximum distance traveled in 2D array.
        Inputs: arr_2D (np.ndarray): 2D array.
                mesh (Mesh): Mesh object containing information about the grid.
        Outputs: tuple: A tuple containing:
                - dis_max (float): Maximum distance traveled by the tracer (furthest coordinate).
                - dis_min (float): Minimum distance traveled by the tracer (closest coordinate).
        """

        # Create an empty list to store max and min indices.
        list_max = []
        list_min = []
        
        # Iterate through each column of the 2D tracer array.
        for i in range(arr_2D.shape[0]):
            # Find indices where the tracer value is not zero and greater than or equal to 0.1.
            list_indi = np.where((arr_2D[i]!=0))[0]
            
            # If there are valid indices found.
            if list_indi.size > 0: 
                # Get the last valid index (furthest position) and append to list.
                last_max_indi = list_indi[-1]
                list_max.append(last_max_indi)
                # Get the first valid index (closest position) and append to list.
                last_min_indi = list_indi[0]
                list_min.append(last_min_indi)
                
        # Find the max and min index and calculate the distance of it.
        try:
            indi_max = np.max(list_max)
            indi_min = np.min(list_min)

            # Cell size is calculated as total domain width divided by number of cell.
            dis_max = (((np.abs(mesh.bounds[2]) + np.abs(mesh.bounds[3])) / arr_2D.shape[1]) * (indi_max + 0.5)) + mesh.bounds[2]
            dis_min = (((np.abs(mesh.bounds[2]) + np.abs(mesh.bounds[3])) / arr_2D.shape[1]) * (indi_min)) + mesh.bounds[2]

        # Handle the case where no valid indices were found (empty lists).
        except ValueError:
            dis_max = 0
            dis_min = 0
                
        return dis_max, dis_min


    def export_tracking_data(self):
        """
        Export tracking data for shock fronts and dense cloud regions to a CSV file.
        Outputs:
            None: Results are saved to a CSV file at the specified output location.
        """
        # Initialize simulation reader for the main simulation data.
        simulation = SimulationReader(self.input_route_vtk_folder)

        # Get time steps
        time_myr = simulation.read_time_steps()[1]

        # Read initial mesh and variables to calculate maximum density.
        mesh_vtk0 = simulation.read_mesh(0, "data")
        rho_3D0 = simulation.read_variable("rho")
        tr1_3D0 = simulation.read_variable("tr1")
        
        # Calculate cloud density
        rho_clo_3D0 = rho_3D0 * tr1_3D0

        # Find maximum density value for threshold calculations.
        max_rho = np.max(rho_clo_3D0)
        
        # Initialize lists to store tracking data.
        list_cloud_dense=[]
        list_mach_max=[]
        list_mach_min=[]

        for i in range(0, int(len(time_myr))):
            
            print(f"Tracking in simulation file {i}.")

            # Initialize readers for simulation and Mach number data.
            sim_object = SimulationReader(self.input_route_vtk_folder)
            mach_object = SimulationReader(self.input_folder_mach)

             # Read mesh data for both simulation and Mach number datasets.
            mesh_sim = sim_object.read_mesh(i, "data")
            mesh_mach = mach_object.read_mesh(i, "HD_shock", vtk_or_vtr="vtr")

            # Read 3D variables.
            tr1_3D = sim_object.read_variable("tr1")
            rho_3D = sim_object.read_variable("rho")
            mach_clo_3D = mach_object.read_variable("mach_clo")

            # Initialize mesh operations for slicing.
            init_oper_sim= MeshOperations(mesh_sim)
            init_oper_mach= MeshOperations(mesh_mach)
            # Convert 3D data to 2D slice.
            tr1_2D = init_oper_sim.slicing(tr1_3D)
            rho_2D = init_oper_sim.slicing(rho_3D)
            mach_clo_2D = init_oper_mach.slicing(mach_clo_3D)
            
            # Apply thresholds to identify regions of interest.
            tr1_2D_cut = np.where(tr1_2D >= (self.threshold_tracer), tr1_2D, 0.0)               # Tracer threshold for cloud material.
            rho_2D_cut = np.where(rho_2D >= (max_rho/self.threshold_rho), 1.0, 0.0)             # Density threshold for dense regions.
            mach_clo_2D_cut = np.where(mach_clo_2D >= self.threshold_mach, mach_clo_2D, 0.0)    # Mach number threshold for internal shocks.

            # Combine density and tracer thresholds to identify dense cloud material.
            rho_clo_2D_cut = rho_2D_cut * tr1_2D_cut
            # Apply cloud mask to Mach number data to identify shocks in cloud material.
            mach_clo_2D_cut = mach_clo_2D_cut * rho_clo_2D_cut

            # Calculate the distances.
            max_dist_tr1 = self.distance_travel_2D(rho_clo_2D_cut, mesh_sim)[0]
            max_dist_mach_clo = self.distance_travel_2D(mach_clo_2D_cut, mesh_mach)[0]
            min_dist_mach_clo = self.distance_travel_2D(mach_clo_2D_cut, mesh_mach)[1]

            # Store results for each time step.
            list_cloud_dense.append(max_dist_tr1)
            list_mach_max.append(max_dist_mach_clo)
            list_mach_min.append(min_dist_mach_clo)

        # Create a DataFrame with the tracking data.
        data = {'last_mach': list_mach_max, 'first_mach': list_mach_min, 'last_cloud_dense': list_cloud_dense}
        df = pd.DataFrame(data)

        # Define output CSV filename.
        csv_filename = self.name_tracking_output + '.csv'

        # Save DataFrame to CSV file.
        df.to_csv(f'{self.output_tracking_dir}/{csv_filename}', index=False)

        print(f'Se ha creado el archivo CSV: {self.output_tracking_dir}/{csv_filename}')

##############################################################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--input_route_vtk_folder', type=str, required=False, help="Directory name of the VTK simulation files.")
    parser.add_argument('--input_folder_mach', type=str, required=False, help="Directory name of the VTR mach files.")
    parser.add_argument('--name_mach_file', type=str, required=False, help="Base name for the VTR mach file.")
    parser.add_argument('--output_tracking_dir', type=str, required=False, help="Directory name of the output tracking CSV file.")
    parser.add_argument('--name_tracking_output', type=str, required=False, help="Name for output tracking CSV file.")
    parser.add_argument('--threshold_tracer', type=float, required=False, help="Threshold number for tracer.")
    parser.add_argument('--threshold_rho', type=float, required=False, help="Threshold number for density.")
    parser.add_argument('--threshold_mach', type=float, required=False, help="Threshold number for mach number.")
    
    args = parser.parse_args()

    # Use configuration from a YAML file if no command-line arguments are provided.
    if args.input_route_vtk_folder==None and args.input_folder_mach==None and args.name_tracking_output==None and args.output_tracking_dir==None and args.name_mach_file==None and args.threshold_tracer==None and args.threshold_rho==None and args.threshold_mach==None:
        with open("config.yaml", "r") as f:
            config = yaml.safe_load(f)
        
        input_route_vtk_folder = config['route_vtk_folder']
        input_folder_mach = config['calculate_shock_cell']['output_vtr_shocks_folder']
        name_mach_file = config['calculate_shock_cell']['prefix_name_file']
        output_tracking_dir = config['tracking_shocks']['output_tracking_dir']
        name_tracking_output = config['tracking_shocks']['name_tracking_output'] + '_' + str(config['tracking_shocks']['threshold_tracer']) + '_' + str(config['tracking_shocks']['threshold_rho']) + '_' + str(config['tracking_shocks']['threshold_mach'])
        threshold_tracer = config['tracking_shocks']['threshold_tracer']
        threshold_rho = float(config['tracking_shocks']['threshold_rho'])
        threshold_mach = float(config['tracking_shocks']['threshold_mach'])

        last_index = last_index_position(input_route_vtk_folder, input_folder_mach, name_tracking_output, output_tracking_dir, name_mach_file, threshold_tracer, threshold_rho, threshold_mach)
        last_index.export_tracking_data()
    
    else:
        # Use command-line arguments if provided.
        last_index = last_index_position(args.input_route_vtk_folder, args.input_folder_mach, args.name_tracking_output, args.output_tracking_dir, args.name_mach_file, args.threshold_tracer, args.threshold_rho, args.threshold_mach)
        last_index.export_tracking_data()