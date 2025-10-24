"""
Kepler Equation Solver and Analytical Orbit Calculations
=========================================================

This module implements analytical solutions for orbital mechanics problems,
including Kepler's equation solver and position/velocity calculations from
orbital elements.

Author: Generated for SMAD Moon Communications Constellation Study
Date: October 2025
"""

import numpy as np
from .Constants import MU_EARTH, KEPLER_TOLERANCE
from .coe_sv import coe_from_sv

# =============================================================================
# KEPLER'S EQUATION SOLVER
# =============================================================================

def solve_kepler_equation(mean_anomaly, eccentricity, tolerance=KEPLER_TOLERANCE, max_iterations=1000):
    """
    Solve Kepler's equation using Newton-Raphson method.
    
    Kepler's equation: M = E - e*sin(E)
    Where: M = mean anomaly, E = eccentric anomaly, e = eccentricity
    
    Args:
        mean_anomaly (float): Mean anomaly [rad]
        eccentricity (float): Orbital eccentricity [-]
        tolerance (float): Convergence tolerance
        max_iterations (int): Maximum number of iterations
        
    Returns:
        float: Eccentric anomaly [rad]
        
    Raises:
        ValueError: If convergence is not achieved
    """
    # Initial guess for eccentric anomaly
    if eccentricity < 0.8:
        E = mean_anomaly  # Good initial guess for low eccentricity
    else:
        # Better initial guess for high eccentricity
        E = np.pi if mean_anomaly > np.pi else mean_anomaly
    
    for i in range(max_iterations):
        # Calculate function value and derivative
        f = E - eccentricity * np.sin(E) - mean_anomaly
        f_prime = 1 - eccentricity * np.cos(E)
        
        # Newton-Raphson update
        E_new = E - f / f_prime
        
        # Check convergence
        if abs(E_new - E) < tolerance:
            return E_new
            
        E = E_new
    
    raise ValueError(f"Kepler equation did not converge after {max_iterations} iterations")

def mean_anomaly_from_time(time, mean_motion, mean_anomaly_0=0.0):
    """
    Calculate mean anomaly at given time.
    
    Args:
        time (float): Time since epoch [s]
        mean_motion (float): Mean motion [rad/s]
        mean_anomaly_0 (float): Initial mean anomaly [rad]
        
    Returns:
        float: Mean anomaly [rad]
    """
    M = mean_anomaly_0 + mean_motion * time
    # Normalize to [0, 2π]
    return M % (2 * np.pi)

def true_anomaly_from_eccentric(eccentric_anomaly, eccentricity):
    """
    Calculate true anomaly from eccentric anomaly.
    
    Args:
        eccentric_anomaly (float): Eccentric anomaly [rad]
        eccentricity (float): Orbital eccentricity [-]
        
    Returns:
        float: True anomaly [rad]
    """
    # Using the relationship: tan(ν/2) = √((1+e)/(1-e)) * tan(E/2)
    sqrt_factor = np.sqrt((1 + eccentricity) / (1 - eccentricity))
    true_anomaly = 2 * np.arctan(sqrt_factor * np.tan(eccentric_anomaly / 2))
    
    # Ensure correct quadrant
    if eccentric_anomaly > np.pi:
        true_anomaly += 2 * np.pi
    
    return true_anomaly % (2 * np.pi)

# =============================================================================
# ORBITAL POSITION AND VELOCITY CALCULATIONS
# =============================================================================

def orbital_radius(semimajor_axis, eccentricity, true_anomaly):
    """
    Calculate orbital radius from true anomaly.
    
    Args:
        semimajor_axis (float): Semi-major axis [km]
        eccentricity (float): Orbital eccentricity [-]
        true_anomaly (float): True anomaly [rad]
        
    Returns:
        float: Orbital radius [km]
    """
    return semimajor_axis * (1 - eccentricity**2) / (1 + eccentricity * np.cos(true_anomaly))

def position_velocity_from_elements(semimajor_axis, eccentricity, inclination, 
                                  raan, arg_periapsis, true_anomaly, mu=MU_EARTH):
    """
    Calculate position and velocity vectors from orbital elements.
    
    Args:
        semimajor_axis (float): Semi-major axis [km]
        eccentricity (float): Orbital eccentricity [-]
        inclination (float): Inclination [rad]
        raan (float): Right ascension of ascending node [rad]
        arg_periapsis (float): Argument of periapsis [rad]
        true_anomaly (float): True anomaly [rad]
        mu (float): Gravitational parameter [km³/s²]
        
    Returns:
        tuple: (position_vector, velocity_vector) in ECI frame [km, km/s]
    """
    # Calculate orbital radius
    r = orbital_radius(semimajor_axis, eccentricity, true_anomaly)
    
    # Specific orbital energy and angular momentum
    specific_energy = -mu / (2 * semimajor_axis)
    h = np.sqrt(mu * semimajor_axis * (1 - eccentricity**2))
    
    # Position and velocity in perifocal frame
    r_pf = np.array([
        r * np.cos(true_anomaly),
        r * np.sin(true_anomaly),
        0.0
    ])
    
    v_pf = np.array([
        -np.sqrt(mu / (semimajor_axis * (1 - eccentricity**2))) * np.sin(true_anomaly),
        np.sqrt(mu / (semimajor_axis * (1 - eccentricity**2))) * (eccentricity + np.cos(true_anomaly)),
        0.0
    ])
    
    # Transformation matrix from perifocal to ECI
    T_pf_to_eci = transformation_matrix_pf_to_eci(inclination, raan, arg_periapsis)
    
    # Transform to ECI frame
    r_eci = T_pf_to_eci @ r_pf
    v_eci = T_pf_to_eci @ v_pf
    
    return r_eci, v_eci

def transformation_matrix_pf_to_eci(inclination, raan, arg_periapsis):
    """
    Create transformation matrix from perifocal to ECI frame.
    
    Args:
        inclination (float): Inclination [rad]
        raan (float): Right ascension of ascending node [rad]
        arg_periapsis (float): Argument of periapsis [rad]
        
    Returns:
        np.ndarray: 3x3 transformation matrix
    """
    # Rotation matrices
    cos_raan = np.cos(raan)
    sin_raan = np.sin(raan)
    cos_inc = np.cos(inclination)
    sin_inc = np.sin(inclination)
    cos_arg = np.cos(arg_periapsis)
    sin_arg = np.sin(arg_periapsis)
    
    # Combined transformation matrix
    T = np.array([
        [cos_raan * cos_arg - sin_raan * sin_arg * cos_inc,
         -cos_raan * sin_arg - sin_raan * cos_arg * cos_inc,
         sin_raan * sin_inc],
        [sin_raan * cos_arg + cos_raan * sin_arg * cos_inc,
         -sin_raan * sin_arg + cos_raan * cos_arg * cos_inc,
         -cos_raan * sin_inc],
        [sin_arg * sin_inc,
         cos_arg * sin_inc,
         cos_inc]
    ])
    
    return T

# =============================================================================
# ANALYTICAL ORBIT PROPAGATION
# =============================================================================

def analytical_propagation(orbital_elements, time_array, mu=MU_EARTH):
    """
    Analytically propagate orbit over time array using Kepler's laws.
    
    This function implements the complete analytical solution to the two-body problem.
    It uses the following sequence:
    1. Extract orbital elements and calculate mean motion
    2. Convert initial true anomaly to initial mean anomaly 
    3. For each time point:
       - Calculate mean anomaly at that time (linear with time)
       - Solve Kepler's equation to find eccentric anomaly
       - Convert eccentric anomaly to true anomaly
       - Calculate position and velocity from true anomaly
    
    Args:
        orbital_elements (dict): Dictionary containing orbital elements
        time_array (np.ndarray): Array of times [s] since initial epoch
        mu (float): Gravitational parameter [km³/s²]
        
    Returns:
        tuple: (positions, velocities) arrays of shape (N, 3) in ECI frame
    """
    # STEP 1: Extract orbital elements from dictionary
    a = orbital_elements['semimajor_axis']     # Semi-major axis [km]
    e = orbital_elements['eccentricity']       # Eccentricity [-]
    i = orbital_elements['inclination']        # Inclination [rad]
    raan = orbital_elements['raan']            # Right ascension of ascending node [rad]
    arg_p = orbital_elements['arg_periapsis']  # Argument of periapsis [rad]
    nu_0 = orbital_elements['true_anomaly_0']  # Initial true anomaly [rad]
    
    # STEP 2: Calculate mean motion using Kepler's 3rd law
    # n = √(μ/a³) - this is how fast the satellite moves on average
    n = np.sqrt(mu / a**3)  # [rad/s]
    
    # STEP 3: Convert initial true anomaly to initial mean anomaly
    # This is necessary because mean anomaly is what varies linearly with time
    
    # First convert true anomaly to eccentric anomaly using the relationship:
    # tan(ν/2) = √((1+e)/(1-e)) * tan(E/2)
    # Solving for E: E = 2*arctan(√((1-e)/(1+e)) * tan(ν/2))
    E_0 = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(nu_0 / 2))
    
    # Handle quadrant issues - if true anomaly > π, eccentric anomaly needs adjustment
    if nu_0 > np.pi:
        E_0 += np.pi
    
    # Convert eccentric anomaly to mean anomaly using Kepler's equation:
    # M = E - e*sin(E)
    M_0 = E_0 - e * np.sin(E_0)  # Initial mean anomaly [rad]
    
    # STEP 4: Initialize arrays to store results
    positions = []   # Will store [x, y, z] positions in ECI frame
    velocities = []  # Will store [vx, vy, vz] velocities in ECI frame
    
    # STEP 5: For each time point, calculate position and velocity
    for t in time_array:
        
        # 5a) Calculate mean anomaly at time t
        # Mean anomaly increases linearly with time: M(t) = M₀ + n*t
        M = mean_anomaly_from_time(t, n, M_0)
        
        # 5b) Solve Kepler's equation: M = E - e*sin(E) for E
        # This is transcendental and requires iterative solution (Newton-Raphson)
        E = solve_kepler_equation(M, e)
        
        # 5c) Convert eccentric anomaly back to true anomaly
        # This tells us the satellite's actual angular position in its orbit
        nu = true_anomaly_from_eccentric(E, e)
        
        # 5d) Calculate 3D position and velocity vectors
        # This uses the true anomaly and orbital elements to get ECI coordinates
        r, v = position_velocity_from_elements(a, e, i, raan, arg_p, nu, mu)
        
        # 5e) Store the results
        positions.append(r)
        velocities.append(v)
    
    # STEP 6: Convert lists to numpy arrays and return
    return np.array(positions), np.array(velocities)

def analytical_true_anomaly_propagation(orbital_elements, time_array, mu=MU_EARTH):
    """
    Analytically propagate true anomaly over time array using Kepler's laws.
    
    This function provides the exact true anomaly as a function of time for 
    comparison with numerical integration results. It follows the same steps
    as analytical_propagation but returns only the true anomaly values.
    
    Args:
        orbital_elements (dict): Dictionary containing orbital elements
        time_array (np.ndarray): Array of times [s] since initial epoch
        mu (float): Gravitational parameter [km³/s²]
        
    Returns:
        np.ndarray: True anomaly values [rad] at each time point
    """
    # Extract orbital elements
    a = orbital_elements['semimajor_axis']
    e = orbital_elements['eccentricity']
    nu_0 = orbital_elements['true_anomaly_0']
    
    # Calculate mean motion
    n = np.sqrt(mu / a**3)
    
    # Convert initial true anomaly to initial mean anomaly
    E_0 = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(nu_0 / 2))
    if nu_0 > np.pi:
        E_0 += np.pi
    M_0 = E_0 - e * np.sin(E_0)
    
    # Calculate true anomaly for each time point
    true_anomalies = []
    
    for t in time_array:
        # Calculate mean anomaly at time t
        M = mean_anomaly_from_time(t, n, M_0)
        
        # Solve Kepler's equation for eccentric anomaly
        E = solve_kepler_equation(M, e)
        
        # Convert to true anomaly
        nu = true_anomaly_from_eccentric(E, e)
        
        true_anomalies.append(nu)
    
    return np.array(true_anomalies)

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def orbital_elements_from_state_vector(r_vec, v_vec, mu=MU_EARTH):
    """
    Calculate orbital elements from state vector using the robust coe_from_sv implementation.
    
    This function is now a wrapper around the more comprehensive coe_from_sv function
    from coe_sv.py, which handles special cases (circular, equatorial orbits) better
    and follows Algorithm 4.1 from orbital mechanics literature.
    
    Args:
        r_vec (np.ndarray): Position vector [km]
        v_vec (np.ndarray): Velocity vector [km/s]
        mu (float): Gravitational parameter [km³/s²]
        
    Returns:
        dict: Dictionary containing orbital elements with consistent naming
    """
    # Use the robust implementation from coe_sv.py
    # Returns: [h, e, RA, incl, w, TA, a]
    coe_list = coe_from_sv(r_vec, v_vec, mu)
    
    # Extract values from the returned list
    h = coe_list[0]           # Specific angular momentum [km²/s]
    e = coe_list[1]           # Eccentricity [-]
    raan = coe_list[2]        # Right ascension of ascending node [rad]
    i = coe_list[3]           # Inclination [rad]
    arg_p = coe_list[4]       # Argument of periapsis [rad]
    nu = coe_list[5]          # True anomaly [rad]
    a = coe_list[6]           # Semi-major axis [km]
    
    # Return dictionary with consistent naming convention used throughout the project
    return {
        'semimajor_axis': a,
        'eccentricity': e,
        'inclination': i,
        'raan': raan,
        'arg_periapsis': arg_p,
        'true_anomaly': nu,
        'specific_angular_momentum': h  # Additional parameter from robust implementation
    }

if __name__ == "__main__":
    # Test the Kepler equation solver
    print("Testing Kepler equation solver...")
    
    # Test case: e=0.9, M=π/4
    e_test = 0.9
    M_test = np.pi / 4
    
    E_solution = solve_kepler_equation(M_test, e_test)
    nu_solution = true_anomaly_from_eccentric(E_solution, e_test)
    
    print(f"Eccentricity: {e_test}")
    print(f"Mean anomaly: {M_test:.6f} rad ({np.degrees(M_test):.2f}°)")
    print(f"Eccentric anomaly: {E_solution:.6f} rad ({np.degrees(E_solution):.2f}°)")
    print(f"True anomaly: {nu_solution:.6f} rad ({np.degrees(nu_solution):.2f}°)")
    
    # Verify solution
    M_check = E_solution - e_test * np.sin(E_solution)
    print(f"Verification: M_calculated = {M_check:.6f} rad (error: {abs(M_check - M_test):.2e})")