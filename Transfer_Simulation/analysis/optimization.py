"""
Multi-parameter optimization functions for lunar transfer missions
"""

import numpy as np
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import DEFAULT_TARGET_PERIGEE_ALTITUDE
from operations.earth_operations import calculate_earth_departure_delta_v
from operations.trajectory_calculations import lunar_trajectory_calculations
from operations.lunar_operations import lunar_soi_calculations, calculate_lunar_soi_transit_time, hyperbolic_to_elliptical_conversion


def multi_parameter_optimization(num_steps=360):
    """
    Comprehensive optimization function that varies R0, V0, gamma0, and lambda1
    to find the optimal combination that minimizes mission delta-V.
    
    Parameters:
    - num_steps: Number of steps for the optimization (minimum 360)
    
    Parameter ranges:
    - R0: Fixed at 1.05 DU
    - V0: 1.372 to 2.0 DU/TU
    - gamma0: 0 to 90 degrees
    - lambda1: 0 to 360 degrees
    
    Returns:
    - Dictionary containing optimization results and optimal parameters
    """
    print("MULTI-PARAMETER OPTIMIZATION STUDY")
    print("=" * 70)
    print(f"Running optimization with {num_steps} total steps...")
    
    # Fixed parameters
    R0 = 1.05  # DU - Fixed as requested
    
    # Variable parameter ranges
    V0_min, V0_max = 1.372, 2.0  # DU/TU
    gamma0_min, gamma0_max = 0, 90  # degrees
    lambda1_min, lambda1_max = 0, 360  # degrees
    
    # Calculate steps per parameter (cube root for 3D parameter space)
    steps_per_param = int(np.ceil(num_steps ** (1/3)))
    total_combinations = steps_per_param ** 3
    
    print(f"Parameter ranges:")
    print(f"  R0: {R0} DU (fixed)")
    print(f"  V0: {V0_min} to {V0_max} DU/TU ({steps_per_param} steps)")
    print(f"  gamma0: {gamma0_min}° to {gamma0_max}° ({steps_per_param} steps)")
    print(f"  lambda1: {lambda1_min}° to {lambda1_max}° ({steps_per_param} steps)")
    print(f"Total combinations: {total_combinations}")
    print()
    
    # Create parameter grids
    V0_values = np.linspace(V0_min, V0_max, steps_per_param)
    gamma0_values = np.linspace(gamma0_min, gamma0_max, steps_per_param)
    lambda1_values = np.linspace(lambda1_min, lambda1_max, steps_per_param)
    
    # Storage for results
    results = []
    best_result = None
    min_delta_v = float('inf')
    
    # Progress tracking
    calculation_count = 0
    successful_calculations = 0
    
    print("Starting optimization calculations...")
    
    for i, V0 in enumerate(V0_values):
        for j, gamma0 in enumerate(gamma0_values):
            for k, lambda1 in enumerate(lambda1_values):
                calculation_count += 1
                
                try:
                    # Calculate Earth departure delta-V
                    departure_results = calculate_earth_departure_delta_v(R0, V0, verbose=False)
                    earth_departure_dv = departure_results['delta_v_departure_kms']
                    
                    # Perform geocentric trajectory calculations
                    geo_results = lunar_trajectory_calculations(R0, V0, gamma0, lambda1, verbose=False)
                    
                    # Perform lunar SOI calculations
                    lunar_results = lunar_soi_calculations(
                        geo_results['r1'], 
                        geo_results['v1'], 
                        geo_results['phi1_deg'],
                        lambda1,
                        geo_results['gamma1_deg'],
                        verbose=False
                    )
                    
                    # Calculate SOI transit time
                    soi_transit_results = calculate_lunar_soi_transit_time(lunar_results, verbose=False)
                    
                    # Calculate hyperbolic to elliptical conversion
                    elliptical_results = hyperbolic_to_elliptical_conversion(
                        lunar_results, 
                        target_perigee_altitude_km=DEFAULT_TARGET_PERIGEE_ALTITUDE, 
                        verbose=False
                    )
                    
                    # Calculate total mission delta-V and time
                    total_delta_v = earth_departure_dv + elliptical_results['total_delta_v']
                    total_time = geo_results['tof_hours'] + soi_transit_results['soi_transit_time_hours']
                    
                    # Store result
                    result = {
                        'R0': R0,
                        'V0': V0,
                        'gamma0': gamma0,
                        'lambda1': lambda1,
                        'earth_departure_dv': earth_departure_dv,
                        'lunar_dv': elliptical_results['total_delta_v'],
                        'total_delta_v': total_delta_v,
                        'total_time_hours': total_time,
                        'geo_tof_hours': geo_results['tof_hours'],
                        'soi_transit_hours': soi_transit_results['soi_transit_time_hours']
                    }
                    
                    results.append(result)
                    successful_calculations += 1
                    
                    # Check if this is the best result so far
                    if total_delta_v < min_delta_v:
                        min_delta_v = total_delta_v
                        best_result = result.copy()
                        
                        print(f"New best found! ΔV = {total_delta_v:.3f} km/s")
                        print(f"  V0 = {V0:.3f} DU/TU, γ0 = {gamma0:.1f}°, λ1 = {lambda1:.1f}°")
                    
                except Exception as e:
                    # Skip failed calculations silently
                    continue
                
                # Progress update every 100 calculations
                if calculation_count % 100 == 0:
                    success_rate = (successful_calculations / calculation_count) * 100
                    print(f"Progress: {calculation_count}/{total_combinations} ({calculation_count/total_combinations*100:.1f}%) - Success rate: {success_rate:.1f}%")
    
    print(f"\nOptimization complete!")
    print(f"Total calculations attempted: {calculation_count}")
    print(f"Successful calculations: {successful_calculations}")
    print(f"Success rate: {(successful_calculations/calculation_count)*100:.1f}%")
    
    if best_result is None:
        print("ERROR: No successful calculations found!")
        return None
    
    # Display optimal results
    print("\n" + "=" * 70)
    print("OPTIMAL MISSION PARAMETERS")
    print("=" * 70)
    print(f"Minimum Mission Delta-V: {best_result['total_delta_v']:.3f} km/s")
    print(f"Optimal Parameters:")
    print(f"  R0 = {best_result['R0']:.3f} DU")
    print(f"  V0 = {best_result['V0']:.3f} DU/TU")
    print(f"  γ0 = {best_result['gamma0']:.1f}°")
    print(f"  λ1 = {best_result['lambda1']:.1f}°")
    print(f"\nMission Breakdown:")
    print(f"  Earth Departure ΔV: {best_result['earth_departure_dv']:.3f} km/s")
    print(f"  Lunar Operations ΔV: {best_result['lunar_dv']:.3f} km/s")
    print(f"  Total Mission ΔV: {best_result['total_delta_v']:.3f} km/s")
    print(f"  Total Mission Time: {best_result['total_time_hours']:.1f} hours")
    print(f"    - Geocentric Transfer: {best_result['geo_tof_hours']:.1f} hours")
    print(f"    - SOI Transit: {best_result['soi_transit_hours']:.1f} hours")
    
    # Analyze parameter sensitivity
    if len(results) > 10:
        print("\n" + "=" * 70)
        print("PARAMETER SENSITIVITY ANALYSIS")
        print("=" * 70)
        
        results_array = np.array(results)
        delta_v_array = np.array([r['total_delta_v'] for r in results])
        
        print(f"Delta-V Statistics:")
        print(f"  Minimum: {np.min(delta_v_array):.3f} km/s")
        print(f"  Maximum: {np.max(delta_v_array):.3f} km/s")
        print(f"  Mean: {np.mean(delta_v_array):.3f} km/s")
        print(f"  Standard Deviation: {np.std(delta_v_array):.3f} km/s")
        print(f"  Range: {np.max(delta_v_array) - np.min(delta_v_array):.3f} km/s")
    
    return {
        'optimal_result': best_result,
        'all_results': results,
        'statistics': {
            'total_calculations': calculation_count,
            'successful_calculations': successful_calculations,
            'success_rate': (successful_calculations/calculation_count)*100,
            'min_delta_v': np.min([r['total_delta_v'] for r in results]),
            'max_delta_v': np.max([r['total_delta_v'] for r in results]),
            'mean_delta_v': np.mean([r['total_delta_v'] for r in results]),
            'std_delta_v': np.std([r['total_delta_v'] for r in results])
        }
    }


def save_optimization_results(results, filename="optimization_results.npz"):
    """
    Save optimization results to a file for later analysis
    
    Parameters:
    - results: Dictionary returned from multi_parameter_optimization
    - filename: Name of file to save results to
    """
    if results is None or 'all_results' not in results:
        print("No valid results to save!")
        return
    
    # Extract data for saving
    all_results = results['all_results']
    
    # Convert to arrays
    R0_values = [r['R0'] for r in all_results]
    V0_values = [r['V0'] for r in all_results]
    gamma0_values = [r['gamma0'] for r in all_results]
    lambda1_values = [r['lambda1'] for r in all_results]
    earth_departure_dv = [r['earth_departure_dv'] for r in all_results]
    lunar_dv = [r['lunar_dv'] for r in all_results]
    total_delta_v = [r['total_delta_v'] for r in all_results]
    total_time_hours = [r['total_time_hours'] for r in all_results]
    
    # Create results directory if it doesn't exist
    results_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
    os.makedirs(results_dir, exist_ok=True)
    
    # Save to compressed numpy file
    filepath = os.path.join(results_dir, filename)
    np.savez_compressed(filepath,
                       R0_values=R0_values,
                       V0_values=V0_values,
                       gamma0_values=gamma0_values,
                       lambda1_values=lambda1_values,
                       earth_departure_dv=earth_departure_dv,
                       lunar_dv=lunar_dv,
                       total_delta_v=total_delta_v,
                       total_time_hours=total_time_hours,
                       optimal_result=results['optimal_result'],
                       statistics=results['statistics'])
    
    print(f"Optimization results saved to: {filepath}")


def load_optimization_results(filename="optimization_results.npz"):
    """
    Load previously saved optimization results
    
    Parameters:
    - filename: Name of file to load results from
    
    Returns:
    - Dictionary containing loaded optimization results
    """
    filepath = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results", filename)
    
    try:
        data = np.load(filepath, allow_pickle=True)
        
        # Reconstruct results dictionary
        results = {
            'optimal_result': data['optimal_result'].item(),
            'statistics': data['statistics'].item(),
            'all_results': []
        }
        
        # Reconstruct individual results
        for i in range(len(data['V0_values'])):
            result = {
                'R0': data['R0_values'][i],
                'V0': data['V0_values'][i],
                'gamma0': data['gamma0_values'][i],
                'lambda1': data['lambda1_values'][i],
                'earth_departure_dv': data['earth_departure_dv'][i],
                'lunar_dv': data['lunar_dv'][i],
                'total_delta_v': data['total_delta_v'][i],
                'total_time_hours': data['total_time_hours'][i]
            }
            results['all_results'].append(result)
        
        print(f"Loaded optimization results from: {filepath}")
        print(f"Total results: {len(results['all_results'])}")
        
        return results
        
    except FileNotFoundError:
        print(f"File not found: {filepath}")
        return None
    except Exception as e:
        print(f"Error loading results: {e}")
        return None