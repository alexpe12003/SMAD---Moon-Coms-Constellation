"""
Analysis and optimization functions for lunar transfer missions
"""

import numpy as np
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import *
from operations.earth_operations import calculate_earth_departure_delta_v
from operations.trajectory_calculations import lunar_trajectory_calculations
from operations.lunar_operations import lunar_soi_calculations, calculate_lunar_soi_transit_time, hyperbolic_to_elliptical_conversion


def parametric_study_lambda1():
    """
    Perform a parametric study of lambda1 from 0 to 360 degrees with 5-degree steps
    and analyze the relationship between lambda1, total mission delta-V, and time of flight
    
    Returns:
    - Tuple of (lambda1_values, total_delta_v_values, time_of_flight_values)
    """
    print("PARAMETRIC STUDY: LAMBDA1 vs MISSION PARAMETERS")
    print("=" * 60)
    
    # Fixed initial conditions
    R0 = DEFAULT_PARKING_ORBIT_RADIUS  # DU - Both parking orbit radius AND transfer trajectory perigee
    V0 = DEFAULT_TRANSFER_VELOCITY  # DU/TU
    gamma0 = DEFAULT_FLIGHT_PATH_ANGLE  # degrees
    
    # Lambda1 range: 0 to 360 degrees with 5-degree steps
    lambda1_values = np.arange(LAMBDA1_RANGE[0], LAMBDA1_RANGE[1] + 1, LAMBDA1_STEP)
    
    # Arrays to store results
    total_delta_v_values = []  # Complete mission delta-V including Earth departure
    time_of_flight_values = []  # Complete mission time including SOI transit
    valid_lambda1_values = []
    
    print(f"Running calculations for {len(lambda1_values)} lambda1 values...")
    print("This may take a moment...")
    
    # Calculate Earth departure delta-V once (it's constant for all lambda1)
    departure_results = calculate_earth_departure_delta_v(R0, V0, verbose=False)
    earth_departure_dv = departure_results['delta_v_departure_kms']
    
    for i, lambda1 in enumerate(lambda1_values):
        try:
            # Perform calculations without verbose output
            geo_results = lunar_trajectory_calculations(R0, V0, gamma0, lambda1, verbose=False)
            lunar_results = lunar_soi_calculations(
                geo_results['r1'], 
                geo_results['v1'], 
                geo_results['phi1_deg'],
                lambda1,
                geo_results['gamma1_deg'],
                verbose=False
            )
            soi_transit_results = calculate_lunar_soi_transit_time(lunar_results, verbose=False)
            elliptical_results = hyperbolic_to_elliptical_conversion(lunar_results, target_perigee_altitude_km=DEFAULT_TARGET_PERIGEE_ALTITUDE, verbose=False)
            
            # Store complete mission results
            complete_delta_v = earth_departure_dv + elliptical_results['total_delta_v']
            complete_time = geo_results['tof_hours'] + soi_transit_results['soi_transit_time_hours']
            
            total_delta_v_values.append(complete_delta_v)
            time_of_flight_values.append(complete_time)
            valid_lambda1_values.append(lambda1)
            
            # Progress indicator
            if (i + 1) % 10 == 0:
                print(f"Completed {i + 1}/{len(lambda1_values)} calculations...")
                
        except Exception as e:
            print(f"Warning: Calculation failed for lambda1 = {lambda1}°: {str(e)}")
            continue
    
    print(f"Successfully calculated {len(valid_lambda1_values)} data points.")
    print(f"Earth departure delta-V (constant): {earth_departure_dv:.3f} km/s")
    
    return valid_lambda1_values, total_delta_v_values, time_of_flight_values


def find_optimal_lambda1(lambda1_values, delta_v_values, tof_values):
    """
    Find and display optimal lambda1 values for different criteria
    
    Parameters:
    - lambda1_values: List of lambda1 values in degrees
    - delta_v_values: List of corresponding delta-V values in km/s
    - tof_values: List of corresponding time of flight values in hours
    
    Returns:
    - Dictionary containing optimal solutions
    """
    print("\n" + "=" * 60)
    print("OPTIMAL LAMBDA1 ANALYSIS")
    print("=" * 60)
    
    # Convert to numpy arrays for easier manipulation
    lambda1_array = np.array(lambda1_values)
    delta_v_array = np.array(delta_v_values)
    time_array = np.array(tof_values)
    
    # Find minimum delta-V
    min_delta_v_idx = np.argmin(delta_v_array)
    min_delta_v = delta_v_array[min_delta_v_idx]
    optimal_lambda1_delta_v = lambda1_array[min_delta_v_idx]
    corresponding_time = time_array[min_delta_v_idx]
    
    print(f"MINIMUM COMPLETE MISSION DELTA-V CRITERION:")
    print(f"  Optimal λ₁ = {optimal_lambda1_delta_v}°")
    print(f"  ΔV = {min_delta_v:.3f} km/s")
    print(f"  Mission Time = {corresponding_time:.1f} hours")
    
    # Find minimum time of flight
    min_time_idx = np.argmin(time_array)
    min_time = time_array[min_time_idx]
    optimal_lambda1_time = lambda1_array[min_time_idx]
    corresponding_delta_v = delta_v_array[min_time_idx]
    
    print(f"\nMINIMUM MISSION TIME CRITERION:")
    print(f"  Optimal λ₁ = {optimal_lambda1_time}°")
    print(f"  Mission Time = {min_time:.1f} hours")
    print(f"  ΔV = {corresponding_delta_v:.3f} km/s")
    
    # Find balanced solution (normalized weighted sum)
    # Normalize both criteria to 0-1 range
    delta_v_normalized = (delta_v_array - np.min(delta_v_array)) / (np.max(delta_v_array) - np.min(delta_v_array))
    time_normalized = (time_array - np.min(time_array)) / (np.max(time_array) - np.min(time_array))
    
    # Equal weighting for balanced solution
    combined_score = 0.5 * delta_v_normalized + 0.5 * time_normalized
    balanced_idx = np.argmin(combined_score)
    balanced_lambda1 = lambda1_array[balanced_idx]
    balanced_delta_v = delta_v_array[balanced_idx]
    balanced_time = time_array[balanced_idx]
    
    print(f"\nBALANCED SOLUTION (Equal weight ΔV and Time):")
    print(f"  Optimal λ₁ = {balanced_lambda1}°")
    print(f"  ΔV = {balanced_delta_v:.3f} km/s")
    print(f"  Mission Time = {balanced_time:.1f} hours")
    
    # Statistics
    print(f"\nSTATISTICS:")
    print(f"  ΔV Range: {np.min(delta_v_array):.3f} - {np.max(delta_v_array):.3f} km/s")
    print(f"  Time Range: {np.min(time_array):.1f} - {np.max(time_array):.1f} hours")
    print(f"  ΔV Standard Deviation: {np.std(delta_v_array):.3f} km/s")
    print(f"  Time Standard Deviation: {np.std(time_array):.1f} hours")
    
    return {
        'min_delta_v': {'lambda1': optimal_lambda1_delta_v, 'delta_v': min_delta_v, 'time': corresponding_time},
        'min_time': {'lambda1': optimal_lambda1_time, 'delta_v': corresponding_delta_v, 'time': min_time},
        'balanced': {'lambda1': balanced_lambda1, 'delta_v': balanced_delta_v, 'time': balanced_time}
    }


def find_optimal_lambda1_complete(lambda1_values, complete_delta_v_values, complete_time_values):
    """
    Find and display optimal lambda1 values for complete mission criteria
    
    Parameters:
    - lambda1_values: List of lambda1 values in degrees
    - complete_delta_v_values: List of complete mission delta-V values in km/s
    - complete_time_values: List of complete mission time values in hours
    
    Returns:
    - Dictionary containing optimal solutions
    """
    print("\n" + "=" * 60)
    print("OPTIMAL LAMBDA1 ANALYSIS (COMPLETE MISSION)")
    print("=" * 60)
    
    # Convert to numpy arrays for easier manipulation
    lambda1_array = np.array(lambda1_values)
    delta_v_array = np.array(complete_delta_v_values)
    time_array = np.array(complete_time_values)
    
    # Find minimum delta-V
    min_delta_v_idx = np.argmin(delta_v_array)
    min_delta_v = delta_v_array[min_delta_v_idx]
    optimal_lambda1_delta_v = lambda1_array[min_delta_v_idx]
    corresponding_time = time_array[min_delta_v_idx]
    
    print(f"MINIMUM COMPLETE MISSION DELTA-V CRITERION:")
    print(f"  Optimal λ₁ = {optimal_lambda1_delta_v}°")
    print(f"  ΔV = {min_delta_v:.3f} km/s")
    print(f"  Mission Time = {corresponding_time:.1f} hours")
    
    # Find minimum time of flight
    min_time_idx = np.argmin(time_array)
    min_time = time_array[min_time_idx]
    optimal_lambda1_time = lambda1_array[min_time_idx]
    corresponding_delta_v = delta_v_array[min_time_idx]
    
    print(f"\nMINIMUM COMPLETE MISSION TIME CRITERION:")
    print(f"  Optimal λ₁ = {optimal_lambda1_time}°")
    print(f"  Mission Time = {min_time:.1f} hours")
    print(f"  ΔV = {corresponding_delta_v:.3f} km/s")
    
    # Find balanced solution
    delta_v_normalized = (delta_v_array - np.min(delta_v_array)) / (np.max(delta_v_array) - np.min(delta_v_array))
    time_normalized = (time_array - np.min(time_array)) / (np.max(time_array) - np.min(time_array))
    
    # Equal weighting for balanced solution
    combined_score = 0.5 * delta_v_normalized + 0.5 * time_normalized
    balanced_idx = np.argmin(combined_score)
    balanced_lambda1 = lambda1_array[balanced_idx]
    balanced_delta_v = delta_v_array[balanced_idx]
    balanced_time = time_array[balanced_idx]
    
    print(f"\nBALANCED SOLUTION (Equal weight ΔV and Time):")
    print(f"  Optimal λ₁ = {balanced_lambda1}°")
    print(f"  ΔV = {balanced_delta_v:.3f} km/s")
    print(f"  Mission Time = {balanced_time:.1f} hours")
    
    # Statistics
    print(f"\nSTATISTICS:")
    print(f"  ΔV Range: {np.min(delta_v_array):.3f} - {np.max(delta_v_array):.3f} km/s")
    print(f"  Time Range: {np.min(time_array):.1f} - {np.max(time_array):.1f} hours")
    print(f"  ΔV Standard Deviation: {np.std(delta_v_array):.3f} km/s")
    print(f"  Time Standard Deviation: {np.std(time_array):.1f} hours")
    
    return {
        'min_delta_v': {'lambda1': optimal_lambda1_delta_v, 'delta_v': min_delta_v, 'time': corresponding_time},
        'min_time': {'lambda1': optimal_lambda1_time, 'delta_v': corresponding_delta_v, 'time': min_time},
        'balanced': {'lambda1': balanced_lambda1, 'delta_v': balanced_delta_v, 'time': balanced_time}
    }