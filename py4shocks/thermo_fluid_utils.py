# thermo_fluid_utils.py
# modified by: Bryan J. Pinargote (bryanpinargote2000@gmail.com)
# mar. 2025

##############################################################################################################

# Import required libraries
import numpy as np

##############################################################################################################

def sound_speed(gamma, prs_3D, rho_3D):
    """
    Calculate the speed of sound in a medium given its pressure and density.
    Inputs: gamma (float): The adiabatic index of the gas.
            prs_3D (numpy.ndarray): A 3D array of pressure values.
            rho_3D (numpy.ndarray): A 3D array of density values.
    Outputs: numpy.ndarray: A 3D array of speed of sound values.
    """
    c_so_3D = np.sqrt(gamma * prs_3D / rho_3D)
    return c_so_3D


def pression_gradient_modulus(prs_3D, dx):
    """
    Calculate the normalized modulus of the pressure gradient.
    This function computes the gradient of the pressure field in three dimensions, then calculates the modulus of this gradient vector. The result is normalized by dividing by the pressure values and multiplying by the grid spacing.
    Inputs: prs_cgs3D (numpy.ndarray): A 3D array of pressure values in CGS units.
            dx (float): The grid spacing in the same units as the pressure.
    Outputs: mdl_prs_gra (numpy.ndarray): A 3D array containing the normalized modulus of the pressure gradient.
    """
    # Calculate the gradient of pressure in all three dimensions
    prs_dx, prs_dy, prs_dz = np.gradient(prs_3D, dx)
    
    # Compute the modulus of the gradient vector
    modulus_gradient_prs = np.sqrt((prs_dx**2) + (prs_dy**2) + (prs_dz**2))
    
    # Normalize the modulus of the gradient by the pressure values and grid spacing
    mdl_prs_gra = (modulus_gradient_prs * dx) / prs_3D
    
    return mdl_prs_gra


def diver_vel(comp1, comp2, comp3, dx, dy, dz):
    """
    Calculate the divergence of a vector field.
    This function computes the divergence of a vector field given its three components and the spatial increments in each direction.
    Inputs: comp1 (numpy.ndarray): Component of the vector field along the x-axis.
            comp2 (numpy.ndarray): Component of the vector field along the y-axis.
            comp3 (numpy.ndarray): Component of the vector field along the z-axis.
            dx (float): Grid spacing in the x-direction.
            dy (float): Grid spacing in the y-direction.
            dz (float): Grid spacing in the z-direction.
    Outputs: tuple: A tuple containing the divergence calculated along each axis and the total one.
    """
    # Calculate the gradient of the first component along the x-axis
    div1 = np.gradient(comp3, dz, axis=0)
    # Calculate the gradient of the second component along the y-axis
    div2 = np.gradient(comp2, dy, axis=1)
    # Calculate the gradient of the third component along the z-axis
    div3 = np.gradient(comp1, dx, axis=2)

    divT = div1 + div2 + div3
    
    return div1, div2, div3, divT


def temperature_func(prs_3D, rho_3D, mu):
    """
    Calculate the temperature from pressure and density arrays using the formula:
    T = (amu * mu * P) / (rho * k_B)
    where mu is the atomic mass unit, m is the mean particle mass, P is the pressure, rho is the density, and k_B is the Boltzmann constant.
    Inputs: prs_cgs3D (numpy.ndarray): 3D array of pressure values in CGS units.
            rho_cgs3D (numpy.ndarray): 3D array of density values in CGS units.
            mu (float): composition - mean particle mass.
    Outputs: temp_K3D (numpy.ndarray): 3D array of temperature values in Kelvin.
    """
    # Define constants
    amu = 1.66053886e-24  # Atomic mass unit in grams
    k_B = 1.3806505e-16   # Boltzmann constant in CGS units (erg/K)
    
    # Calculate temperature in Kelvin
    temp_K3D = ((amu * mu * prs_3D) / (rho_3D * k_B))

    return temp_K3D
