"""
Main application for lunar transfer trajectory analysis

This is the main entry point for the lunar transfer analysis program.
It provides a clean interface for running either single calculations or
parametric studies of lunar transfer trajectories.
"""

# Import all required modules
import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from core.interface import display_analysis_menu, get_user_input, display_organized_mission_summary, get_optimization_parameters
from operations.earth_operations import calculate_earth_departure_delta_v
from operations.trajectory_calculations import lunar_trajectory_calculations
from operations.lunar_operations import lunar_soi_calculations, calculate_lunar_soi_transit_time, hyperbolic_to_elliptical_conversion
from analysis.analysis import parametric_study_lambda1, find_optimal_lambda1
from analysis.plotting import create_parametric_plots
from analysis.optimization import multi_parameter_optimization



def perform_single_calculation():
    """
    Perform a single lunar transfer calculation with user-defined lambda1
    
    Returns:
    - Tuple containing all calculation results
    """
    # Get user inputs
    user_inputs = get_user_input()
    
    # Extract variables
    R0 = user_inputs['R0']
    V0 = user_inputs['V0']
    gamma0 = user_inputs['gamma0']
    lambda1 = user_inputs['lambda1']
    
    # Perform complete mission analysis with detailed output
    departure_results = calculate_earth_departure_delta_v(R0, V0, verbose=True)
    geo_results = lunar_trajectory_calculations(R0, V0, gamma0, lambda1, verbose=True)
    lunar_results = lunar_soi_calculations(
        geo_results['r1'], 
        geo_results['v1'], 
        geo_results['phi1_deg'],
        lambda1,
        geo_results['gamma1_deg'],
        verbose=True
    )
    soi_transit_results = calculate_lunar_soi_transit_time(lunar_results, verbose=True)
    circular_results = hyperbolic_to_elliptical_conversion(
        lunar_results, 
        verbose=True
    )
    
    # Display complete organized mission summary
    display_organized_mission_summary(
        user_inputs, departure_results, geo_results, 
        lunar_results, soi_transit_results, circular_results
    )
    
    return user_inputs, departure_results, geo_results, lunar_results, soi_transit_results, circular_results


def perform_parametric_study():
    """
    Perform a parametric study of lambda1 from 0-360 degrees
    
    Returns:
    - Tuple containing parametric study results
    """
    # Run parametric study
    lambda1_values, delta_v_values, tof_values = parametric_study_lambda1()
    
    # Create plots
    create_parametric_plots(lambda1_values, delta_v_values, tof_values)
    
    # Find optimal values
    optimal_results = find_optimal_lambda1(lambda1_values, delta_v_values, tof_values)
    
    return lambda1_values, delta_v_values, tof_values, optimal_results


def perform_multi_parameter_optimization():
    """
    Perform multi-parameter optimization of V0, gamma0, and lambda1
    
    Returns:
    - Dictionary containing optimization results
    """
    # Get optimization parameters from user
    opt_params = get_optimization_parameters()
    
    # Run optimization
    print("\nStarting multi-parameter optimization...")
    results = multi_parameter_optimization(num_steps=opt_params['num_steps'])
    
    # Note: File saving capability has been removed as requested
    if results is None:
        print("Optimization failed.")
    
    return results


def main():
    """
    Main function to run the complete lunar transfer trajectory analysis
    """
    # Display menu and get user choice
    choice = display_analysis_menu()
    
    if choice == '1':
        # Single calculation
        results = perform_single_calculation()
        return results
    elif choice == '2':
        # Parametric study
        results = perform_parametric_study()
        return results
    else:
        # Multi-parameter optimization
        results = perform_multi_parameter_optimization()
        return results


if __name__ == "__main__":
    results = main()