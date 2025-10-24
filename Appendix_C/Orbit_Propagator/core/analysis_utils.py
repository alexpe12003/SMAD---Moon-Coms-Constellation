"""
Analysis and Plotting Utilities for Orbit Propagation Studies
============================================================

This module provides advanced analysis tools and plotting utilities
for studying orbital mechanics simulations, including convergence
analysis, error metrics, and visualization tools.

Author: Generated for SMAD Moon Communications Constellation Study
Date: October 2025
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from scipy import interpolate
from .Constants import RESULTS_DIR
import os

# =============================================================================
# ERROR ANALYSIS FUNCTIONS
# =============================================================================

def calculate_position_error_metrics(numerical_pos, analytical_pos, times):
    """
    Calculate comprehensive position error metrics.
    
    Args:
        numerical_pos (np.ndarray): Numerical solution positions [N, 3]
        analytical_pos (np.ndarray): Analytical solution positions [N, 3]
        times (np.ndarray): Time array [N]
        
    Returns:
        dict: Dictionary containing error metrics
    """
    # Calculate position errors
    position_errors = np.linalg.norm(numerical_pos - analytical_pos, axis=1)
    
    # Error metrics
    metrics = {
        'max_error': np.max(position_errors),
        'mean_error': np.mean(position_errors),
        'rms_error': np.sqrt(np.mean(position_errors**2)),
        'final_error': position_errors[-1],
        'error_growth_rate': (position_errors[-1] - position_errors[0]) / (times[-1] - times[0]),
        'relative_max_error': np.max(position_errors) / np.mean(np.linalg.norm(analytical_pos, axis=1)),
        'error_time_series': position_errors
    }
    
    return metrics

def calculate_velocity_error_metrics(numerical_vel, analytical_vel, times):
    """
    Calculate comprehensive velocity error metrics.
    
    Args:
        numerical_vel (np.ndarray): Numerical solution velocities [N, 3]
        analytical_vel (np.ndarray): Analytical solution velocities [N, 3]
        times (np.ndarray): Time array [N]
        
    Returns:
        dict: Dictionary containing error metrics
    """
    # Calculate velocity errors
    velocity_errors = np.linalg.norm(numerical_vel - analytical_vel, axis=1)
    
    # Error metrics
    metrics = {
        'max_error': np.max(velocity_errors),
        'mean_error': np.mean(velocity_errors),
        'rms_error': np.sqrt(np.mean(velocity_errors**2)),
        'final_error': velocity_errors[-1],
        'error_growth_rate': (velocity_errors[-1] - velocity_errors[0]) / (times[-1] - times[0]),
        'relative_max_error': np.max(velocity_errors) / np.mean(np.linalg.norm(analytical_vel, axis=1)),
        'error_time_series': velocity_errors
    }
    
    return metrics

def analyze_conservation_properties(positions, velocities, times, mu=398600.4418):
    """
    Analyze conservation of energy and angular momentum.
    
    Args:
        positions (np.ndarray): Position array [N, 3]
        velocities (np.ndarray): Velocity array [N, 3]
        times (np.ndarray): Time array [N]
        mu (float): Gravitational parameter
        
    Returns:
        dict: Conservation analysis results
    """
    # Calculate energy
    r_mag = np.linalg.norm(positions, axis=1)
    v_mag = np.linalg.norm(velocities, axis=1)
    energy = v_mag**2 / 2 - mu / r_mag
    
    # Calculate angular momentum
    h_vec = np.cross(positions, velocities)
    h_mag = np.linalg.norm(h_vec, axis=1)
    
    # Energy conservation analysis
    energy_error = energy - energy[0]
    energy_metrics = {
        'initial_energy': energy[0],
        'final_energy': energy[-1],
        'max_energy_error': np.max(np.abs(energy_error)),
        'relative_energy_error': np.max(np.abs(energy_error)) / np.abs(energy[0]),
        'energy_drift_rate': (energy[-1] - energy[0]) / (times[-1] - times[0]),
        'energy_time_series': energy,
        'energy_error_time_series': energy_error
    }
    
    # Angular momentum conservation analysis
    h_error = h_mag - h_mag[0]
    h_metrics = {
        'initial_h': h_mag[0],
        'final_h': h_mag[-1],
        'max_h_error': np.max(np.abs(h_error)),
        'relative_h_error': np.max(np.abs(h_error)) / h_mag[0],
        'h_drift_rate': (h_mag[-1] - h_mag[0]) / (times[-1] - times[0]),
        'h_time_series': h_mag,
        'h_error_time_series': h_error
    }
    
    return {
        'energy': energy_metrics,
        'angular_momentum': h_metrics
    }

# =============================================================================
# CONVERGENCE ANALYSIS
# =============================================================================

def estimate_convergence_order(step_sizes, errors):
    """
    Estimate convergence order from error vs step size data.
    
    Args:
        step_sizes (np.ndarray): Array of step sizes
        errors (np.ndarray): Corresponding errors
        
    Returns:
        tuple: (convergence_order, r_squared, fit_params)
    """
    # Filter out zero or negative values
    valid_mask = (step_sizes > 0) & (errors > 0)
    
    if np.sum(valid_mask) < 2:
        return None, 0.0, None
    
    log_h = np.log(step_sizes[valid_mask])
    log_err = np.log(errors[valid_mask])
    
    # Linear fit in log-log space: log(error) = p * log(h) + c
    fit_params = np.polyfit(log_h, log_err, 1)
    convergence_order = fit_params[0]
    
    # Calculate R-squared
    log_err_pred = np.polyval(fit_params, log_h)
    ss_res = np.sum((log_err - log_err_pred)**2)
    ss_tot = np.sum((log_err - np.mean(log_err))**2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0.0
    
    return convergence_order, r_squared, fit_params

def richardson_extrapolation(step_sizes, solutions, target_step=None):
    """
    Apply Richardson extrapolation to improve solution accuracy.
    
    Args:
        step_sizes (list): List of step sizes (should be in ratio 2:1)
        solutions (list): Corresponding numerical solutions
        target_step (float): Target step size for extrapolation
        
    Returns:
        np.ndarray: Extrapolated solution
    """
    if len(step_sizes) != 2 or len(solutions) != 2:
        raise ValueError("Richardson extrapolation requires exactly 2 solutions")
    
    h1, h2 = step_sizes
    sol1, sol2 = solutions
    
    # Assume 4th order method (RK4)
    p = 4
    
    # Richardson extrapolation formula
    r = h1 / h2
    extrapolated = sol2 + (sol2 - sol1) / (r**p - 1)
    
    return extrapolated

# =============================================================================
# ADVANCED PLOTTING FUNCTIONS
# =============================================================================

def plot_orbital_elements_evolution(times, positions, velocities, orbital_params, 
                                  title_suffix="", save_path=None):
    """
    Plot evolution of orbital elements over time.
    
    Uses the robust orbital_elements_from_state_vector function which now
    integrates the comprehensive coe_from_sv implementation for better
    handling of special cases (circular, equatorial orbits).
    
    Args:
        times (np.ndarray): Time array
        positions (np.ndarray): Position array [N, 3]
        velocities (np.ndarray): Velocity array [N, 3]
        orbital_params (dict): Reference orbital parameters
        title_suffix (str): Suffix for plot titles
        save_path (str): Path to save plot
    """
    from Kepler import orbital_elements_from_state_vector
    
    # Calculate orbital elements at each time step
    elements_history = []
    for i in range(len(times)):
        elements = orbital_elements_from_state_vector(positions[i], velocities[i])
        elements_history.append(elements)
    
    # Convert to arrays
    a_history = np.array([elem['semimajor_axis'] for elem in elements_history])
    e_history = np.array([elem['eccentricity'] for elem in elements_history])
    i_history = np.array([elem['inclination'] for elem in elements_history])
    
    # Normalize time to orbital periods
    period = orbital_params['orbital_period']
    time_periods = times / period
    
    # Create plots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Semi-major axis
    axes[0, 0].plot(time_periods, a_history, 'b-', linewidth=2)
    axes[0, 0].axhline(orbital_params['semimajor_axis'], color='r', linestyle='--', 
                       label='Theoretical')
    axes[0, 0].set_xlabel('Time (orbital periods)')
    axes[0, 0].set_ylabel('Semi-major axis (km)')
    axes[0, 0].set_title(f'Semi-major Axis Evolution{title_suffix}')
    axes[0, 0].grid(True)
    axes[0, 0].legend()
    
    # Eccentricity
    axes[0, 1].plot(time_periods, e_history, 'g-', linewidth=2)
    axes[0, 1].axhline(orbital_params['eccentricity'], color='r', linestyle='--',
                       label='Theoretical')
    axes[0, 1].set_xlabel('Time (orbital periods)')
    axes[0, 1].set_ylabel('Eccentricity')
    axes[0, 1].set_title(f'Eccentricity Evolution{title_suffix}')
    axes[0, 1].grid(True)
    axes[0, 1].legend()
    
    # Inclination
    axes[1, 0].plot(time_periods, np.degrees(i_history), 'm-', linewidth=2)
    axes[1, 0].axhline(np.degrees(orbital_params['inclination']), color='r', 
                       linestyle='--', label='Theoretical')
    axes[1, 0].set_xlabel('Time (orbital periods)')
    axes[1, 0].set_ylabel('Inclination (degrees)')
    axes[1, 0].set_title(f'Inclination Evolution{title_suffix}')
    axes[1, 0].grid(True)
    axes[1, 0].legend()
    
    # Orbital radius vs true anomaly (orbit shape check)
    r_mag = np.linalg.norm(positions, axis=1)
    nu_history = np.array([elem['true_anomaly'] for elem in elements_history])
    
    # Theoretical orbit
    nu_theory = np.linspace(0, 2*np.pi, 1000)
    r_theory = (orbital_params['semimajor_axis'] * (1 - orbital_params['eccentricity']**2) /
                (1 + orbital_params['eccentricity'] * np.cos(nu_theory)))
    
    axes[1, 1].plot(np.degrees(nu_theory), r_theory, 'r--', label='Theoretical', linewidth=2)
    axes[1, 1].plot(np.degrees(nu_history), r_mag, 'b-', label='Numerical', alpha=0.7)
    axes[1, 1].set_xlabel('True Anomaly (degrees)')
    axes[1, 1].set_ylabel('Orbital Radius (km)')
    axes[1, 1].set_title(f'Orbit Shape Verification{title_suffix}')
    axes[1, 1].grid(True)
    axes[1, 1].legend()
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    plt.show()

def plot_error_analysis_comprehensive(convergence_results, orbital_params, save_path=None):
    """
    Create comprehensive error analysis plots.
    
    Args:
        convergence_results (dict): Results from convergence study
        orbital_params (dict): Orbital parameters
        save_path (str): Path to save plots
    """
    fig = plt.figure(figsize=(16, 12))
    
    # Extract data
    step_sizes = list(convergence_results.keys())
    step_sizes.sort()
    
    max_pos_errors = [convergence_results[h]['max_position_error'] for h in step_sizes]
    max_vel_errors = [convergence_results[h]['max_velocity_error'] for h in step_sizes]
    max_energy_errors = [convergence_results[h]['max_energy_error'] for h in step_sizes]
    comp_times = [convergence_results[h]['computation_time'] for h in step_sizes]
    num_steps = [convergence_results[h]['num_steps'] for h in step_sizes]
    
    # 1. Convergence plot (log-log)
    ax1 = plt.subplot(2, 3, 1)
    plt.loglog(step_sizes, max_pos_errors, 'bo-', linewidth=2, markersize=8, label='Position')
    plt.loglog(step_sizes, max_vel_errors, 'ro-', linewidth=2, markersize=8, label='Velocity')
    
    # Add theoretical convergence lines
    h_ref = np.array(step_sizes)
    error_ref = max_pos_errors[0] * (h_ref / h_ref[0])**4
    plt.loglog(h_ref, error_ref, 'k--', alpha=0.7, label='4th Order')
    
    plt.xlabel('Step Size (s)')
    plt.ylabel('Maximum Error')
    plt.title('Convergence Analysis')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # 2. Computational efficiency
    ax2 = plt.subplot(2, 3, 2)
    plt.loglog(step_sizes, comp_times, 'go-', linewidth=2, markersize=8)
    plt.xlabel('Step Size (s)')
    plt.ylabel('Computation Time (s)')
    plt.title('Computational Cost')
    plt.grid(True, alpha=0.3)
    
    # 3. Accuracy vs efficiency trade-off
    ax3 = plt.subplot(2, 3, 3)
    plt.loglog(max_pos_errors, comp_times, 'mo-', linewidth=2, markersize=8)
    for i, h in enumerate(step_sizes):
        plt.annotate(f'{h}s', (max_pos_errors[i], comp_times[i]), 
                    xytext=(5, 5), textcoords='offset points', fontsize=8)
    plt.xlabel('Maximum Position Error (km)')
    plt.ylabel('Computation Time (s)')
    plt.title('Accuracy vs Efficiency')
    plt.grid(True, alpha=0.3)
    
    # 4. Error evolution for different step sizes
    ax4 = plt.subplot(2, 3, 4)
    period = orbital_params['orbital_period']
    
    for i, h in enumerate(step_sizes[::2]):  # Every other step size
        result = convergence_results[h]
        times = result['times'] / period
        pos_errors = result['position_errors']
        plt.semilogy(times, pos_errors, linewidth=2, label=f'Δt={h}s')
    
    plt.xlabel('Time (orbital periods)')
    plt.ylabel('Position Error (km)')
    plt.title('Error Evolution Over Time')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # 5. Energy conservation
    ax5 = plt.subplot(2, 3, 5)
    for i, h in enumerate(step_sizes[::2]):
        result = convergence_results[h]
        times = result['times'] / period
        energy_errors = result['energy_errors']
        plt.semilogy(times, np.abs(energy_errors), linewidth=2, label=f'Δt={h}s')
    
    plt.xlabel('Time (orbital periods)')
    plt.ylabel('|Energy Error| (km²/s²)')
    plt.title('Energy Conservation')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # 6. Step size efficiency metric
    ax6 = plt.subplot(2, 3, 6)
    efficiency = np.array(comp_times) / np.array(max_pos_errors)
    plt.semilogx(step_sizes, efficiency, 'co-', linewidth=2, markersize=8)
    plt.xlabel('Step Size (s)')
    plt.ylabel('Time/Error Ratio (s/km)')
    plt.title('Efficiency Metric')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    plt.show()

# =============================================================================
# DATA EXPORT FUNCTIONS
# =============================================================================

def export_results_to_csv(convergence_results, filename=None):
    """
    Export convergence results to CSV format.
    
    Args:
        convergence_results (dict): Results from convergence study
        filename (str): Output filename
    """
    if filename is None:
        filename = f"{RESULTS_DIR}/convergence_results.csv"
    
    # Prepare data
    data = []
    for step_size, result in convergence_results.items():
        data.append({
            'step_size_s': step_size,
            'num_steps': result['num_steps'],
            'computation_time_s': result['computation_time'],
            'max_position_error_km': result['max_position_error'],
            'max_velocity_error_km_s': result['max_velocity_error'],
            'max_energy_error_km2_s2': result['max_energy_error'],
            'max_h_error_km2_s': result['max_h_error'],
            'final_position_error_km': result['final_position_error'],
            'final_velocity_error_km_s': result['final_velocity_error']
        })
    
    # Create DataFrame and save
    df = pd.DataFrame(data)
    df.to_csv(filename, index=False)
    
    print(f"Results exported to {filename}")
    return df

def create_analysis_report(simulation_results, orbital_params, filename=None):
    """
    Create a comprehensive analysis report.
    
    Args:
        simulation_results (dict): Complete simulation results
        orbital_params (dict): Orbital parameters
        filename (str): Output filename
    """
    if filename is None:
        filename = f"{RESULTS_DIR}/analysis_report.md"
    
    convergence_data = simulation_results['convergence']
    
    # Estimate convergence order
    step_sizes = np.array(list(convergence_data.keys()))
    max_pos_errors = np.array([convergence_data[h]['max_position_error'] 
                              for h in step_sizes])
    
    conv_order, r_squared, _ = estimate_convergence_order(step_sizes, max_pos_errors)
    
    with open(filename, 'w') as f:
        f.write("# Orbit B Numerical Integration Analysis Report\n\n")
        
        f.write("## Orbital Parameters\n")
        for key, value in orbital_params.items():
            f.write(f"- **{key.replace('_', ' ').title()}**: {value}\n")
        
        f.write(f"\n## Convergence Analysis\n")
        f.write(f"- **Estimated Convergence Order**: {conv_order:.2f}\n")
        f.write(f"- **R-squared**: {r_squared:.4f}\n")
        f.write(f"- **Expected for RK4**: 4.0\n\n")
        
        f.write("## Step Size Analysis\n")
        f.write("| Step Size (s) | Steps | Time (s) | Max Pos Error (km) | Max Vel Error (km/s) |\n")
        f.write("|---------------|-------|----------|--------------------|-----------------------|\n")
        
        for h in sorted(step_sizes):
            result = convergence_data[h]
            f.write(f"| {h:>11.1f} | {result['num_steps']:>5d} | "
                   f"{result['computation_time']:>8.3f} | "
                   f"{result['max_position_error']:>18.3e} | "
                   f"{result['max_velocity_error']:>21.3e} |\n")
        
        f.write(f"\n## Recommendations\n")
        
        # Find optimal step size based on error/time trade-off
        efficiencies = []
        for h in step_sizes:
            result = convergence_data[h]
            efficiency = result['computation_time'] / result['max_position_error']
            efficiencies.append((h, efficiency))
        
        optimal_h = min(efficiencies, key=lambda x: x[1])[0]
        f.write(f"- **Optimal Step Size**: {optimal_h} s (best accuracy/time ratio)\n")
        
        # Performance recommendations
        fastest_h = min(step_sizes)
        most_accurate = max(step_sizes)
        
        f.write(f"- **For Speed**: Use {fastest_h} s step size\n")
        f.write(f"- **For Accuracy**: Use {most_accurate} s step size\n")
        f.write(f"- **Balanced**: Use {optimal_h} s step size\n")
    
    print(f"Analysis report saved to {filename}")

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def ensure_results_directory():
    """Ensure results directory exists."""
    os.makedirs(RESULTS_DIR, exist_ok=True)

if __name__ == "__main__":
    # Test analysis functions
    print("Analysis utilities module loaded successfully!")
    print("Available functions:")
    print("- calculate_position_error_metrics")
    print("- calculate_velocity_error_metrics") 
    print("- analyze_conservation_properties")
    print("- estimate_convergence_order")
    print("- plot_orbital_elements_evolution")
    print("- plot_error_analysis_comprehensive")
    print("- export_results_to_csv")
    print("- create_analysis_report")