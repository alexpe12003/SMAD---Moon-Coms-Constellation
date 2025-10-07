"""
Main application for lunar transfer trajectory analysis

This is the main entry point for the lunar transfer analysis program.
It provides a clean interface for running either single calculations or
parametric studies of lunar transfer trajectories.
"""

# Import all required modules
from interface import display_analysis_menu, get_user_input, display_organized_mission_summary
from earth_operations import calculate_earth_departure_delta_v
from trajectory_calculations import lunar_trajectory_calculations
from lunar_operations import lunar_soi_calculations, calculate_lunar_soi_transit_time, hyperbolic_to_elliptical_conversion
from analysis import parametric_study_lambda1, find_optimal_lambda1
from plotting import create_parametric_plots
from config import DEFAULT_TARGET_PERIGEE_ALTITUDE


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
    elliptical_results = hyperbolic_to_elliptical_conversion(
        lunar_results, 
        target_perigee_altitude_km=DEFAULT_TARGET_PERIGEE_ALTITUDE, 
        verbose=True
    )
    
    # Display complete organized mission summary
    display_organized_mission_summary(
        user_inputs, departure_results, geo_results, 
        lunar_results, soi_transit_results, elliptical_results
    )
    
    return user_inputs, departure_results, geo_results, lunar_results, soi_transit_results, elliptical_results


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
    else:
        # Parametric study
        results = perform_parametric_study()
        return results


if __name__ == "__main__":
    results = main()