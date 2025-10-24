#!/usr/bin/env python3
"""
Fixed True Anomaly Calculation Functions
"""

import numpy as np

def true_anomaly_from_eccentric_old(eccentric_anomaly, eccentricity):
    """Original function for comparison."""
    sqrt_factor = np.sqrt((1 + eccentricity) / (1 - eccentricity))
    true_anomaly = 2 * np.arctan(sqrt_factor * np.tan(eccentric_anomaly / 2))
    
    if eccentric_anomaly > np.pi:
        true_anomaly += 2 * np.pi
    
    return true_anomaly % (2 * np.pi)

def true_anomaly_from_eccentric_fixed(eccentric_anomaly, eccentricity):
    """
    FIXED: Calculate true anomaly from eccentric anomaly.
    
    The issue with the original function was improper quadrant handling.
    This version properly handles all quadrants and ensures continuity.
    
    Args:
        eccentric_anomaly (float): Eccentric anomaly [rad]
        eccentricity (float): Orbital eccentricity [-]
        
    Returns:
        float: True anomaly [rad] in range [0, 2π]
    """
    # Use the robust relationship: cos(ν) = (cos(E) - e) / (1 - e*cos(E))
    # This avoids the tan(ν/2) formula which has quadrant issues
    
    cos_E = np.cos(eccentric_anomaly)
    cos_nu = (cos_E - eccentricity) / (1 - eccentricity * cos_E)
    
    # Clamp to avoid numerical issues
    cos_nu = np.clip(cos_nu, -1.0, 1.0)
    
    # Calculate the principal value (0 to π)
    nu = np.arccos(cos_nu)
    
    # Determine correct quadrant:
    # If E is in [0, π], then ν is in [0, π]
    # If E is in [π, 2π], then ν is in [π, 2π]
    if eccentric_anomaly > np.pi:
        nu = 2*np.pi - nu
    
    return nu

def calculate_true_anomaly_from_state_fixed(position, velocity, mu):
    """
    FIXED: Calculate true anomaly from position and velocity vectors.
    
    The issue with the original function was in the quadrant determination logic.
    This version uses a more robust approach based on the flight path angle.
    
    Args:
        position (np.ndarray): Position vector(s) [km]
        velocity (np.ndarray): Velocity vector(s) [km/s]  
        mu (float): Gravitational parameter [km³/s²]
        
    Returns:
        float or np.ndarray: True anomaly [rad]
    """
    # Handle both single vectors and arrays
    single_vector = position.ndim == 1
    if single_vector:
        position = position.reshape(1, 3)
        velocity = velocity.reshape(1, 3)
    
    # Calculate orbital elements to get true anomaly directly
    # This is more robust than trying to reconstruct it from eccentricity vector
    
    results = []
    for i in range(len(position)):
        r_vec = position[i]
        v_vec = velocity[i]
        
        # Calculate specific angular momentum vector
        h_vec = np.cross(r_vec, v_vec)
        h = np.linalg.norm(h_vec)
        
        # Calculate orbital radius and speed
        r = np.linalg.norm(r_vec)
        v = np.linalg.norm(v_vec)
        
        # Calculate specific energy and eccentricity
        energy = v**2/2 - mu/r
        a = -mu/(2*energy)  # Semi-major axis
        
        # Calculate eccentricity vector and magnitude
        e_vec = ((v**2 - mu/r) * r_vec - np.dot(r_vec, v_vec) * v_vec) / mu
        e = np.linalg.norm(e_vec)
        
        # Calculate true anomaly using the dot product method
        if e < 1e-8:  # Nearly circular orbit
            # For circular orbits, use angular position from reference direction
            # Use the direction of the eccentricity vector as reference (even if small)
            if np.linalg.norm(e_vec) < 1e-10:
                # Pure circular - use x-axis as reference
                cos_nu = r_vec[0] / r
                sin_nu = r_vec[1] / r  # Assuming motion in x-y plane
            else:
                cos_nu = np.dot(e_vec, r_vec) / (np.linalg.norm(e_vec) * r)
                # Use velocity direction to determine quadrant
                h_cross_r = np.cross(h_vec, r_vec)
                sin_nu = np.dot(e_vec, h_cross_r) / (np.linalg.norm(e_vec) * np.linalg.norm(h_cross_r))
        else:
            # For eccentric orbits, use standard formula
            cos_nu = np.dot(e_vec, r_vec) / (e * r)
            # Determine quadrant using radial velocity
            radial_velocity = np.dot(r_vec, v_vec) / r
            sin_nu = np.sqrt(1 - cos_nu**2)
            if radial_velocity < 0:  # Moving toward periapsis
                sin_nu = -sin_nu
        
        # Calculate true anomaly using atan2 for proper quadrant
        nu = np.arctan2(sin_nu, cos_nu)
        
        # Ensure result is in [0, 2π]
        if nu < 0:
            nu += 2*np.pi
            
        results.append(nu)
    
    return results[0] if single_vector else np.array(results)

if __name__ == "__main__":
    # Test the fixed functions
    print("Testing FIXED true anomaly functions:")
    print("=" * 60)
    
    # Test 1: Eccentric to true anomaly conversion
    print("1. Testing eccentric to true anomaly conversion:")
    e = 0.9
    E_values = np.linspace(0, 2*np.pi, 9)
    
    print(f"{'E (deg)':>8} {'ν old (deg)':>12} {'ν fixed (deg)':>14} {'Δν fixed':>12}")
    print("-" * 60)
    
    nu_prev = None
    for E in E_values:
        nu_old = true_anomaly_from_eccentric_old(E, e) 
        nu_fixed = true_anomaly_from_eccentric_fixed(E, e)
        
        delta_nu = 0 if nu_prev is None else np.degrees(nu_fixed - nu_prev)
        if delta_nu > 180:
            delta_nu -= 360
        elif delta_nu < -180:
            delta_nu += 360
            
        print(f"{np.degrees(E):8.1f} {np.degrees(nu_old):12.1f} {np.degrees(nu_fixed):14.1f} {delta_nu:12.1f}")
        nu_prev = nu_fixed