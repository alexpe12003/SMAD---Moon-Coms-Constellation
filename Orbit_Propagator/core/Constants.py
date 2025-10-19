"""
Orbital Mechanics Constants and Configuration Parameters
========================================================

This module contains all physical constants and configurable parameters
used in the orbit propagation study. Any parameters that need to be
modified for study purposes should be defined here.

Author: Generated for SMAD Moon Communications Constellation Study
Date: October 2025
"""

import numpy as np

# =============================================================================
# PHYSICAL CONSTANTS
# =============================================================================

# Earth gravitational parameter [km³/s²]
MU_EARTH = 398600.4418  # Standard gravitational parameter for Earth

# Earth radius [km]
R_EARTH = 6371.0  # Mean Earth radius

# =============================================================================
# ORBIT B PARAMETERS (Primary Study Orbit)
# =============================================================================

# Orbit B characteristics as specified
ORBIT_B_SEMIMAJOR_AXIS = 70000.0  # km - Semi-major axis
ORBIT_B_ECCENTRICITY = 0.9         # - Eccentricity

# Additional orbital elements for Orbit B (can be modified for study)
ORBIT_B_INCLINATION = 0.0          # rad - Inclination
ORBIT_B_RAAN = 0.0                 # rad - Right Ascension of Ascending Node
ORBIT_B_ARG_PERIAPSIS = 0.0        # rad - Argument of periapsis
ORBIT_B_TRUE_ANOMALY_0 = 0.0       # rad - Initial true anomaly

# =============================================================================
# NUMERICAL INTEGRATION PARAMETERS
# =============================================================================

# Time integration parameters
DEFAULT_TIME_STEP = 100.0            # s - Default integration time step
MIN_TIME_STEP = 1.0                 # s - Minimum allowed time step
MAX_TIME_STEP = 300.0               # s - Maximum allowed time step

# Integration study parameters
INTEGRATION_STEPS = [1.0, 5.0, 10.0, 30.0, 60.0, 100.0, 300.0]  # s - Steps for convergence study
NUM_ORBITS_STUDY = 5                # - Number of orbits to integrate for study

# Convergence tolerances
KEPLER_TOLERANCE = 1e-12            # - Tolerance for Kepler equation solver
POSITION_TOLERANCE = 1e-6           # km - Position convergence tolerance
VELOCITY_TOLERANCE = 1e-9           # km/s - Velocity convergence tolerance

# =============================================================================
# ANALYSIS PARAMETERS
# =============================================================================

# Output parameters
SAVE_RESULTS = True                 # - Whether to save results to files
PLOT_RESULTS = True                 # - Whether to generate plots
RESULTS_DIR = "results"             # - Directory for saving results

# Plotting parameters
FIGURE_SIZE = (12, 8)               # - Default figure size for plots
DPI = 300                           # - Plot resolution
PLOT_STYLE = 'seaborn-v0_8'        # - Matplotlib style

# =============================================================================
# DERIVED PARAMETERS (Calculated from orbit parameters)
# =============================================================================

def get_orbit_b_parameters():
    """
    Calculate derived orbital parameters for Orbit B.
    
    Returns:
        dict: Dictionary containing all orbital parameters
    """
    # Calculate periapsis and apoapsis distances
    r_periapsis = ORBIT_B_SEMIMAJOR_AXIS * (1 - ORBIT_B_ECCENTRICITY)
    r_apoapsis = ORBIT_B_SEMIMAJOR_AXIS * (1 + ORBIT_B_ECCENTRICITY)
    
    # Calculate orbital period
    orbital_period = 2 * np.pi * np.sqrt(ORBIT_B_SEMIMAJOR_AXIS**3 / MU_EARTH)
    
    # Calculate mean motion
    mean_motion = 2 * np.pi / orbital_period
    
    # Calculate specific orbital energy
    specific_energy = -MU_EARTH / (2 * ORBIT_B_SEMIMAJOR_AXIS)
    
    return {
        'semimajor_axis': ORBIT_B_SEMIMAJOR_AXIS,
        'eccentricity': ORBIT_B_ECCENTRICITY,
        'inclination': ORBIT_B_INCLINATION,
        'raan': ORBIT_B_RAAN,
        'arg_periapsis': ORBIT_B_ARG_PERIAPSIS,
        'true_anomaly_0': ORBIT_B_TRUE_ANOMALY_0,
        'r_periapsis': r_periapsis,
        'r_apoapsis': r_apoapsis,
        'orbital_period': orbital_period,
        'mean_motion': mean_motion,
        'specific_energy': specific_energy
    }

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def degrees_to_radians(degrees):
    """Convert degrees to radians."""
    return degrees * np.pi / 180.0

def radians_to_degrees(radians):
    """Convert radians to degrees."""
    return radians * 180.0 / np.pi

def print_orbit_summary():
    """Print a summary of the orbital parameters."""
    params = get_orbit_b_parameters()
    
    print("=" * 60)
    print("ORBIT B PARAMETER SUMMARY")
    print("=" * 60)
    print(f"Semi-major axis:     {params['semimajor_axis']:>12.3f} km")
    print(f"Eccentricity:        {params['eccentricity']:>12.6f}")
    print(f"Periapsis distance:  {params['r_periapsis']:>12.3f} km")
    print(f"Apoapsis distance:   {params['r_apoapsis']:>12.3f} km")
    print(f"Orbital period:      {params['orbital_period']:>12.3f} s ({params['orbital_period']/3600:.3f} hours)")
    print(f"Mean motion:         {params['mean_motion']:>12.9f} rad/s")
    print(f"Specific energy:     {params['specific_energy']:>12.3f} km²/s²")
    print("=" * 60)

if __name__ == "__main__":
    # Print orbit summary when run as main module
    print_orbit_summary()