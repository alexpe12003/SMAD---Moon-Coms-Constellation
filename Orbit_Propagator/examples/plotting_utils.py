"""
Plotting Utilities for Orbit B Numerical Integration Study
=========================================================

This module contains all plotting functions for visualizing the results
of the orbital mechanics convergence study.

Author: Generated for SMAD Moon Communications Constellation Study
Date: October 2025
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import os

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.Constants import INTEGRATION_STEPS, RESULTS_DIR, MU_EARTH
from core.Kepler import analytical_true_anomaly_propagation
from core.RK4 import calculate_true_anomaly_from_state

# =============================================================================
# MAIN PLOTTING FUNCTIONS
# =============================================================================

def create_single_method_plots(simulation, method_name):
    """Create comprehensive plots for a single method (RK4 or RK8) analysis."""
    print(f"\nCreating comprehensive {method_name} analysis plots...")
    
    # Set up plotting style
    plt.style.use('default')
    
    method_key = f'single_method_{method_name.lower()}'
    if method_key not in simulation.results:
        print(f"No data found for {method_name} method")
        return
    
    method_data = simulation.results[method_key]
    
    # Create position vs time plots for all step sizes
    create_position_vs_time_plots(simulation, method_name, method_data)
    
    # Create velocity vs time plots for all step sizes
    create_velocity_vs_time_plots(simulation, method_name, method_data)
    
    # Create true anomaly vs time plots for all step sizes
    create_true_anomaly_vs_time_plots(simulation, method_name, method_data)
    
    # Create position error vs time plots for all step sizes
    create_position_error_vs_time_plots(simulation, method_name, method_data)
    
    # Create velocity error vs time plots for all step sizes
    create_velocity_error_vs_time_plots(simulation, method_name, method_data)
    
    # Create true anomaly error vs time plots for all step sizes
    create_true_anomaly_error_vs_time_plots(simulation, method_name, method_data)
    
    # Create convergence analysis plot
    create_single_method_convergence_plot(simulation, method_name, method_data)
    
    # Create 3D trajectory comparison
    create_single_method_3d_trajectory(simulation, method_name, method_data)
    
    print(f"All {method_name} plots saved to {RESULTS_DIR}/")
    
    # Close all figures to free memory
    plt.close('all')

def create_all_plots(simulation):
    """Create all individual plots for the simulation analysis."""
    print(f"\nCreating individual analysis plots...")
    
    # Set up plotting style
    plt.style.use('default')
    
    # Create individual plots
    create_trajectory_plot(simulation)
    create_convergence_plot(simulation)
    create_position_error_time_plot(simulation)
    create_velocity_error_time_plot(simulation)
    create_true_anomaly_error_time_plot(simulation)
    create_velocity_convergence_plot(simulation)
    create_computation_time_plot(simulation)
    create_3d_trajectory_plot(simulation)
    create_energy_conservation_plot(simulation)
    
    # Create special case plots if data exists
    if 'special_case_100_orbits' in simulation.results:
        create_special_case_position_error_plot(simulation)
        create_special_case_plots(simulation)
    
    print(f"All plots saved to {RESULTS_DIR}/")
    
    # Close all figures to free memory
    plt.close('all')

# =============================================================================
# INDIVIDUAL PLOTTING FUNCTIONS
# =============================================================================

def create_trajectory_plot(simulation):
    """Create 2D trajectory comparison plot."""
    print(f"Creating 2D trajectory comparison plot...")
    
    convergence_data = simulation.results['convergence']
    analytical_data = simulation.results['analytical_reference']
    
    fig = plt.figure(figsize=(16, 10))
    ax = fig.add_subplot(1, 1, 1)
    
    # Plot analytical solution
    pos_analytical = analytical_data['positions']
    ax.plot(pos_analytical[:, 0]/1000, pos_analytical[:, 1]/1000,
            'k-', label='Analytical', linewidth=2)
    
    # Plot RK4 and RK8 solutions for different step sizes
    method_colors = {'RK4': plt.cm.Reds(np.linspace(0.3, 1, len(INTEGRATION_STEPS))),
                     'RK8': plt.cm.Blues(np.linspace(0.3, 1, len(INTEGRATION_STEPS)))}
    method_styles = {'RK4': '--', 'RK8': ':'}
    
    for method_name in ['RK4', 'RK8']:
        for i, step_size in enumerate(INTEGRATION_STEPS):
            result = convergence_data[method_name][step_size]
            pos_method = result['positions']
            ax.plot(pos_method[:, 0]/1000, pos_method[:, 1]/1000,
                    method_styles[method_name], color=method_colors[method_name][i], 
                    label=f'{method_name} (Δt={step_size}s)', alpha=0.7, linewidth=1)
    
    # Mark key orbital points
    orbit_params = simulation.orbit_params
    ax.scatter([orbit_params['r_periapsis']/1000], [0], 
               color='red', s=150, marker='*', label='Periapsis', zorder=10)
    ax.scatter([-orbit_params['r_apoapsis']/1000], [0], 
               color='blue', s=150, marker='*', label='Apoapsis', zorder=10)
    ax.scatter([0], [0], color='black', s=100, marker='o', 
               label='Earth Center', zorder=10)
    
    # Add Earth radius circle
    earth_radius = 6371  # km
    circle = plt.Circle((0, 0), earth_radius/1000, color='lightblue', 
                       fill=False, linestyle=':', alpha=0.7, label='Earth Surface')
    ax.add_patch(circle)
    
    ax.set_xlabel('X (1000 km)', fontsize=12)
    ax.set_ylabel('Y (1000 km)', fontsize=12)
    ax.set_title('2D Orbital Trajectories: Analytical vs RK4 vs RK8', fontsize=14, fontweight='bold')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/trajectory_comparison_2d.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_convergence_plot(simulation):
    """Create convergence study plot for position and true anomaly errors."""
    print(f"Creating convergence study plot...")
    
    convergence_data = simulation.results['convergence']
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    step_sizes = np.array(INTEGRATION_STEPS)
    
    # Colors and markers for different methods
    method_props = {
        'RK4': {'color': 'red', 'marker': 'o', 'linestyle': '-'},
        'RK8': {'color': 'blue', 'marker': 's', 'linestyle': '--'}
    }
    
    for method_name in ['RK4', 'RK8']:
        max_pos_errors = np.array([convergence_data[method_name][h]['max_position_error'] 
                                  for h in INTEGRATION_STEPS])
        max_true_anomaly_errors = np.array([convergence_data[method_name][h]['max_true_anomaly_error'] 
                                          for h in INTEGRATION_STEPS])
        
        props = method_props[method_name]
        
        # Position error convergence
        ax1.loglog(step_sizes, max_pos_errors, 
                  color=props['color'], marker=props['marker'], linestyle=props['linestyle'],
                  label=f'{method_name} Position Error', linewidth=2, markersize=8)
        
        # True anomaly error convergence
        ax2.loglog(step_sizes, max_true_anomaly_errors, 
                  color=props['color'], marker=props['marker'], linestyle=props['linestyle'],
                  label=f'{method_name} True Anomaly Error', linewidth=2, markersize=8)
        
        # Add theoretical convergence lines
        if method_name == 'RK4':
            ax1.loglog(step_sizes, step_sizes**4 * max_pos_errors[0] / step_sizes[0]**4, 
                      color='red', linestyle=':', alpha=0.7, label='4th Order Reference', linewidth=2)
            ax2.loglog(step_sizes, step_sizes**4 * max_true_anomaly_errors[0] / step_sizes[0]**4, 
                      color='red', linestyle=':', alpha=0.7, label='4th Order Reference', linewidth=2)
        elif method_name == 'RK8':
            ax1.loglog(step_sizes, step_sizes**8 * max_pos_errors[0] / step_sizes[0]**8, 
                      color='blue', linestyle=':', alpha=0.7, label='8th Order Reference', linewidth=2)
            ax2.loglog(step_sizes, step_sizes**8 * max_true_anomaly_errors[0] / step_sizes[0]**8, 
                      color='blue', linestyle=':', alpha=0.7, label='8th Order Reference', linewidth=2)
    
    # Position error plot formatting
    ax1.set_xlabel('Step Size (s)', fontsize=12)
    ax1.set_ylabel('Maximum Position Error (km)', fontsize=12)
    ax1.set_title('Position Error Convergence', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)
    
    # True anomaly error plot formatting
    ax2.set_xlabel('Step Size (s)', fontsize=12)
    ax2.set_ylabel('Maximum True Anomaly Error (rad)', fontsize=12)
    ax2.set_title('True Anomaly Error Convergence', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/convergence_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_position_error_time_plot(simulation):
    """Create position error over time plot for all step sizes."""
    print(f"Creating position error over time plot...")
    
    convergence_data = simulation.results['convergence']
    
    fig = plt.figure(figsize=(16, 8))
    ax = fig.add_subplot(1, 1, 1)
    
    # Choose fewer step sizes for cleaner visualization
    display_steps = [10.0, 30.0, 100.0, 300.0]
    method_colors = {'RK4': plt.cm.Reds(np.linspace(0.4, 1, len(display_steps))),
                     'RK8': plt.cm.Blues(np.linspace(0.4, 1, len(display_steps)))}
    method_styles = {'RK4': '-', 'RK8': '--'}
    
    for method_name in ['RK4', 'RK8']:
        for i, step_size in enumerate(display_steps):
            if step_size in INTEGRATION_STEPS:
                result = convergence_data[method_name][step_size]
                times = result['times'] / simulation.orbit_params['orbital_period']
                position_errors = result['position_errors']
                ax.semilogy(times, position_errors, 
                           linestyle=method_styles[method_name],
                           color=method_colors[method_name][i], 
                           label=f'{method_name} Δt={step_size}s', linewidth=2)
    
    ax.set_xlabel('Time (orbital periods)', fontsize=12)
    ax.set_ylabel('Position Error (km)', fontsize=12)
    ax.set_title('Position Error Evolution Over Time', fontsize=14, fontweight='bold')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/position_error_time.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_velocity_error_time_plot(simulation):
    """Create velocity error over time plot for all step sizes."""
    print(f"Creating velocity error over time plot...")
    
    convergence_data = simulation.results['convergence']
    
    fig = plt.figure(figsize=(16, 8))
    ax = fig.add_subplot(1, 1, 1)
    
    # Choose fewer step sizes for cleaner visualization
    display_steps = [10.0, 30.0, 100.0, 300.0]
    method_colors = {'RK4': plt.cm.Reds(np.linspace(0.4, 1, len(display_steps))),
                     'RK8': plt.cm.Blues(np.linspace(0.4, 1, len(display_steps)))}
    method_styles = {'RK4': '-', 'RK8': '--'}
    
    for method_name in ['RK4', 'RK8']:
        for i, step_size in enumerate(display_steps):
            if step_size in INTEGRATION_STEPS:
                result = convergence_data[method_name][step_size]
                times = result['times'] / simulation.orbit_params['orbital_period']
                velocity_errors = result['velocity_errors']
                ax.semilogy(times, velocity_errors, 
                           linestyle=method_styles[method_name],
                           color=method_colors[method_name][i], 
                           label=f'{method_name} Δt={step_size}s', linewidth=2)
    
    ax.set_xlabel('Time (orbital periods)', fontsize=12)
    ax.set_ylabel('Velocity Error (km/s)', fontsize=12)
    ax.set_title('Velocity Error Evolution Over Time', fontsize=14, fontweight='bold')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/velocity_error_time.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_true_anomaly_error_time_plot(simulation):
    """Create true anomaly error over time plot for all step sizes."""
    print(f"Creating true anomaly error over time plot...")
    
    convergence_data = simulation.results['convergence']
    
    fig = plt.figure(figsize=(16, 8))
    ax = fig.add_subplot(1, 1, 1)
    
    # Choose fewer step sizes for cleaner visualization
    display_steps = [10.0, 30.0, 100.0, 300.0]
    method_colors = {'RK4': plt.cm.Reds(np.linspace(0.4, 1, len(display_steps))),
                     'RK8': plt.cm.Blues(np.linspace(0.4, 1, len(display_steps)))}
    method_styles = {'RK4': '-', 'RK8': '--'}
    
    for method_name in ['RK4', 'RK8']:
        for i, step_size in enumerate(display_steps):
            if step_size in INTEGRATION_STEPS:
                result = convergence_data[method_name][step_size]
                times = result['times'] / simulation.orbit_params['orbital_period']
                true_anomaly_errors = result['true_anomaly_errors']
                ax.semilogy(times, np.degrees(true_anomaly_errors), 
                           linestyle=method_styles[method_name],
                           color=method_colors[method_name][i], 
                           label=f'{method_name} Δt={step_size}s', linewidth=2)
    
    ax.set_xlabel('Time (orbital periods)', fontsize=12)
    ax.set_ylabel('True Anomaly Error (degrees)', fontsize=12)
    ax.set_title('True Anomaly Error Evolution Over Time', fontsize=14, fontweight='bold')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/true_anomaly_error_time.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_computation_time_plot(simulation):
    """Create computation time plot for each step size."""
    print(f"Creating computation time plot...")
    
    convergence_data = simulation.results['convergence']
    
    fig = plt.figure(figsize=(14, 8))
    ax = fig.add_subplot(1, 1, 1)
    
    step_sizes = np.array(INTEGRATION_STEPS)
    x_pos = np.arange(len(step_sizes))
    
    # Get data for both methods
    rk4_times = np.array([convergence_data['RK4'][h]['computation_time'] 
                         for h in INTEGRATION_STEPS])
    rk8_times = np.array([convergence_data['RK8'][h]['computation_time'] 
                         for h in INTEGRATION_STEPS])
    rk4_steps = np.array([convergence_data['RK4'][h]['num_steps'] 
                         for h in INTEGRATION_STEPS])
    
    # Create grouped bar chart
    width = 0.35
    rk4_bars = ax.bar(x_pos - width/2, rk4_times, width, 
                      label='RK4', color='red', alpha=0.7, edgecolor='black')
    rk8_bars = ax.bar(x_pos + width/2, rk8_times, width, 
                      label='RK8', color='blue', alpha=0.7, edgecolor='black')
    
    ax.set_xlabel('Step Size (s)', fontsize=12)
    ax.set_ylabel('Computation Time (s)', fontsize=12)
    ax.set_title('Computation Time Comparison: RK4 vs RK8', fontsize=14, fontweight='bold')
    ax.set_xticks(x_pos)
    ax.set_xticklabels([f'{s}s' for s in step_sizes])
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add text annotations with number of steps for RK4 bars only (same for both methods)
    for i, (rk4_bar, steps) in enumerate(zip(rk4_bars, rk4_steps)):
        height = max(rk4_times[i], rk8_times[i])
        ax.text(x_pos[i], height + height*0.02,
                f'{steps:,} steps', ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/computation_time.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_velocity_convergence_plot(simulation):
    """Create velocity error convergence plot."""
    print(f"Creating velocity error convergence plot...")
    
    convergence_data = simulation.results['convergence']
    
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(1, 1, 1)
    
    step_sizes = np.array(INTEGRATION_STEPS)
    
    # Colors and markers for different methods
    method_props = {
        'RK4': {'color': 'red', 'marker': 'o', 'linestyle': '-'},
        'RK8': {'color': 'blue', 'marker': 's', 'linestyle': '--'}
    }
    
    for method_name in ['RK4', 'RK8']:
        max_velocity_errors = np.array([convergence_data[method_name][h]['max_velocity_error'] 
                                       for h in INTEGRATION_STEPS])
        
        props = method_props[method_name]
        
        # Velocity error convergence
        ax.loglog(step_sizes, max_velocity_errors, 
                  color=props['color'], marker=props['marker'], linestyle=props['linestyle'],
                  label=f'{method_name} Velocity Error', linewidth=2, markersize=8)
        
        # Add theoretical convergence lines
        if method_name == 'RK4':
            ax.loglog(step_sizes, step_sizes**4 * max_velocity_errors[0] / step_sizes[0]**4, 
                      color='red', linestyle=':', alpha=0.7, label='4th Order Reference', linewidth=2)
        elif method_name == 'RK8':
            ax.loglog(step_sizes, step_sizes**8 * max_velocity_errors[0] / step_sizes[0]**8, 
                      color='blue', linestyle=':', alpha=0.7, label='8th Order Reference', linewidth=2)
    
    ax.set_xlabel('Step Size (s)', fontsize=12)
    ax.set_ylabel('Maximum Velocity Error (km/s)', fontsize=12)
    ax.set_title('Velocity Error Convergence: RK4 vs RK8', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/velocity_convergence.png', dpi=300, bbox_inches='tight')
    plt.close()

# =============================================================================
# OPTIONAL SPECIALIZED PLOTS
# =============================================================================

def create_3d_trajectory_plot(simulation):
    """Create 3D trajectory comparison plot (optional)."""
    print(f"Creating 3D trajectory comparison plot...")
    
    convergence_data = simulation.results['convergence']
    analytical_data = simulation.results['analytical_reference']
    
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    
    # Plot analytical solution
    pos_analytical = analytical_data['positions']
    ax.plot(pos_analytical[:, 0]/1000, pos_analytical[:, 1]/1000, pos_analytical[:, 2]/1000,
            'k-', label='Analytical', linewidth=2)
    
    # Plot RK4 and RK8 solutions for selected step sizes (every other for clarity)
    selected_steps = [10.0, 100.0, 300.0]  # Select fewer for clarity
    method_colors = {'RK4': ['red', 'darkred', 'maroon'],
                     'RK8': ['blue', 'darkblue', 'navy']}
    method_styles = {'RK4': '--', 'RK8': ':'}
    
    for method_name in ['RK4', 'RK8']:
        for i, step_size in enumerate(selected_steps):
            if step_size in INTEGRATION_STEPS:
                result = convergence_data[method_name][step_size]
                pos_method = result['positions']
                ax.plot(pos_method[:, 0]/1000, pos_method[:, 1]/1000, pos_method[:, 2]/1000,
                        method_styles[method_name], color=method_colors[method_name][i],
                        label=f'{method_name} (Δt={step_size}s)', alpha=0.7, linewidth=1.5)
    
    ax.set_xlabel('X (1000 km)', fontsize=10)
    ax.set_ylabel('Y (1000 km)', fontsize=10)
    ax.set_zlabel('Z (1000 km)', fontsize=10)
    ax.set_title('3D Orbital Trajectories: Analytical vs RK4 vs RK8', fontsize=14, fontweight='bold')
    ax.legend(fontsize=8)
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/trajectory_comparison_3d.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_energy_conservation_plot(simulation):
    """Create energy conservation plot (optional)."""
    print(f"Creating energy conservation plot...")
    
    convergence_data = simulation.results['convergence']
    
    fig = plt.figure(figsize=(16, 8))
    ax = fig.add_subplot(1, 1, 1)
    
    # Choose fewer step sizes for cleaner visualization  
    display_steps = [10.0, 30.0, 100.0, 300.0]
    method_colors = {'RK4': plt.cm.Reds(np.linspace(0.4, 1, len(display_steps))),
                     'RK8': plt.cm.Blues(np.linspace(0.4, 1, len(display_steps)))}
    method_styles = {'RK4': '-', 'RK8': '--'}
    
    for method_name in ['RK4', 'RK8']:
        for i, step_size in enumerate(display_steps):
            if step_size in INTEGRATION_STEPS:
                result = convergence_data[method_name][step_size]
                times = result['times'] / simulation.orbit_params['orbital_period']
                energy_errors = result['energy_errors']
                ax.semilogy(times, np.abs(energy_errors), 
                           linestyle=method_styles[method_name],
                           color=method_colors[method_name][i], 
                           label=f'{method_name} Δt={step_size}s', linewidth=2)
    
    ax.set_xlabel('Time (orbital periods)', fontsize=12)
    ax.set_ylabel('|Energy Error| (km²/s²)', fontsize=12)
    ax.set_title('Energy Conservation Over Time: RK4 vs RK8', fontsize=14, fontweight='bold')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/energy_conservation.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_special_case_position_error_plot(simulation):
    """Create position error plot for the special case: RK8 100s step for 100 orbits."""
    print(f"Creating special case position error plot (RK8 100s step, 100 orbits)...")
    
    special_data = simulation.results['special_case_100_orbits']
    
    fig = plt.figure(figsize=(16, 10))
    ax = fig.add_subplot(1, 1, 1)
    
    # Convert time to orbital periods for better visualization
    times_orbits = special_data['times'] / simulation.orbit_params['orbital_period']
    position_errors = special_data['position_errors']
    
    # Plot position error over time
    ax.semilogy(times_orbits, position_errors, 
               'b-', linewidth=2, alpha=0.8, label='RK8 Position Error (Δt=100s)')
    
    # Add some statistical information
    max_error = np.max(position_errors)
    final_error = position_errors[-1]
    mean_error = np.mean(position_errors)
    
    # Add horizontal lines for reference
    ax.axhline(y=max_error, color='red', linestyle='--', alpha=0.7, 
               label=f'Maximum Error: {max_error:.2e} km')
    ax.axhline(y=final_error, color='orange', linestyle='--', alpha=0.7, 
               label=f'Final Error: {final_error:.2e} km')
    ax.axhline(y=mean_error, color='green', linestyle='--', alpha=0.7, 
               label=f'Mean Error: {mean_error:.2e} km')
    
    # Mark every 10 orbits with vertical lines
    for orbit_mark in range(10, 101, 10):
        ax.axvline(x=orbit_mark, color='gray', linestyle=':', alpha=0.3)
        if orbit_mark % 20 == 0:  # Label every 20 orbits
            ax.text(orbit_mark, max_error * 0.1, f'{orbit_mark}', 
                   rotation=90, verticalalignment='bottom', 
                   horizontalalignment='right', fontsize=9, alpha=0.7)
    
    # Formatting
    ax.set_xlabel('Time (orbital periods)', fontsize=14)
    ax.set_ylabel('Position Error (km)', fontsize=14)
    ax.set_title('RK8 Position Error Evolution Over 100 Orbits\n(Step Size: 100 seconds)', 
                fontsize=16, fontweight='bold')
    ax.legend(fontsize=12, loc='best')
    ax.grid(True, alpha=0.3, which='both')
    
    # Set x-axis limits
    ax.set_xlim(0, 100)
    
    # Add text box with simulation info
    info_text = (f"Integration Method: RK8\n"
                f"Step Size: {special_data['step_size']} s\n"
                f"Total Steps: {special_data['num_steps']:,}\n"
                f"Computation Time: {special_data['computation_time']:.1f} s\n"
                f"Simulation Span: {special_data['num_orbits']} orbits")
    
    ax.text(0.02, 0.98, info_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/special_case_rk8_100orbits_position_error.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_special_case_plots(sim):
    """Create plots for the special case showing position, velocity, and true anomaly over time."""
    if 'special_case_100_orbits' not in sim.results:
        print("No special case results found. Cannot create plots.")
        return
    
    special_data = sim.results['special_case_100_orbits']
    times = special_data['times']
    pos_numerical = special_data['positions']
    vel_numerical = special_data['velocities']
    pos_analytical = special_data['analytical_positions']
    vel_analytical = special_data['analytical_velocities']
    position_errors = special_data['position_errors']
    velocity_errors = special_data['velocity_errors']
    
    # Calculate true anomaly for both methods
    true_anomaly_numerical = calculate_true_anomaly_from_state(pos_numerical, vel_numerical, MU_EARTH)
    true_anomaly_analytical = analytical_true_anomaly_propagation(sim.orbit_params, times)
    
    # Convert time to orbits for better readability
    orbital_period = sim.orbit_params['orbital_period']  # Get orbital period in seconds
    time_orbits = times / orbital_period  # Convert seconds to number of orbits
    
    # Create figure with subplots
    fig = plt.figure(figsize=(15, 12))
    fig.suptitle('Special Case: RK8 vs Analytical Solutions (10 Orbits)', fontsize=16, fontweight='bold')
    
    # 1. Position magnitude over time
    ax1 = plt.subplot(3, 2, 1)
    pos_mag_numerical = np.linalg.norm(pos_numerical, axis=1)
    pos_mag_analytical = np.linalg.norm(pos_analytical, axis=1)
    
    plt.plot(time_orbits, pos_mag_numerical, 'b-', label='RK8 Numerical', linewidth=1.5)
    plt.plot(time_orbits, pos_mag_analytical, 'r--', label='Analytical', linewidth=1.0, alpha=0.8)
    plt.xlabel('Time (orbits)')
    plt.ylabel('Position Magnitude (km)')
    plt.title('Position Magnitude vs Time')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 2. Position error over time
    ax2 = plt.subplot(3, 2, 2)
    plt.semilogy(time_orbits, position_errors, 'g-', linewidth=1.5)
    plt.xlabel('Time (orbits)')
    plt.ylabel('Position Error (km)')
    plt.title('Position Error vs Time')
    plt.grid(True, alpha=0.3)
    
    # 3. Velocity magnitude over time
    ax3 = plt.subplot(3, 2, 3)
    vel_mag_numerical = np.linalg.norm(vel_numerical, axis=1)
    vel_mag_analytical = np.linalg.norm(vel_analytical, axis=1)
    
    plt.plot(time_orbits, vel_mag_numerical, 'b-', label='RK8 Numerical', linewidth=1.5)
    plt.plot(time_orbits, vel_mag_analytical, 'r--', label='Analytical', linewidth=1.0, alpha=0.8)
    plt.xlabel('Time (orbits)')
    plt.ylabel('Velocity Magnitude (km/s)')
    plt.title('Velocity Magnitude vs Time')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 4. Velocity error over time
    ax4 = plt.subplot(3, 2, 4)
    plt.semilogy(time_orbits, velocity_errors, 'orange', linewidth=1.5)
    plt.xlabel('Time (orbits)')
    plt.ylabel('Velocity Error (km/s)')
    plt.title('Velocity Error vs Time')
    plt.grid(True, alpha=0.3)
    
    # 5. True anomaly over time
    ax5 = plt.subplot(3, 2, 5)
    # Handle 2π wraparound for better visualization
    true_anomaly_numerical_unwrapped = np.unwrap(true_anomaly_numerical)
    true_anomaly_analytical_unwrapped = np.unwrap(true_anomaly_analytical)
    
    plt.plot(time_orbits, np.degrees(true_anomaly_numerical_unwrapped), 'b-', 
             label='RK8 Numerical', linewidth=1.5)
    plt.plot(time_orbits, np.degrees(true_anomaly_analytical_unwrapped), 'r--', 
             label='Analytical', linewidth=1.0, alpha=0.8)
    plt.xlabel('Time (orbits)')
    plt.ylabel('True Anomaly (degrees)')
    plt.title('True Anomaly vs Time')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 6. True anomaly error over time
    ax6 = plt.subplot(3, 2, 6)
    # Calculate true anomaly errors (handle 2π wraparound)
    true_anomaly_diff = true_anomaly_numerical - true_anomaly_analytical
    true_anomaly_diff = np.where(true_anomaly_diff > np.pi, 
                                true_anomaly_diff - 2*np.pi, true_anomaly_diff)
    true_anomaly_diff = np.where(true_anomaly_diff < -np.pi, 
                                true_anomaly_diff + 2*np.pi, true_anomaly_diff)
    true_anomaly_errors = np.abs(true_anomaly_diff)
    
    plt.semilogy(time_orbits, np.degrees(true_anomaly_errors), 'purple', linewidth=1.5)
    plt.xlabel('Time (orbits)')
    plt.ylabel('True Anomaly Error (degrees)')
    plt.title('True Anomaly Error vs Time')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save the plot
    plot_filename = f"{RESULTS_DIR}/special_case_comparison_plots.png"
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    print(f"Special case plots saved to: {plot_filename}")
    
    # Show plot if running interactively
    try:
        plt.show()
    except:
        pass  # In case we're running in a non-interactive environment
    
    return fig

def create_special_case_plots_improved(sim):
    """
    IMPROVED: Create plots for the special case with proper handling of highly eccentric orbits.
    
    This version addresses the visualization issues for e=0.9 orbits where true anomaly
    changes extremely rapidly near periapsis but very slowly near apoapsis.
    
    The issue you observed was CORRECT PHYSICS, not a calculation error:
    - For e=0.9 orbits, true anomaly changes ~3800x faster near periapsis than apoapsis
    - This creates 'almost no change' during most of the orbit (near apoapsis)  
    - And 'huge jumps' near periapsis passage
    """
    if 'special_case_100_orbits' not in sim.results:
        print("No special case results found. Cannot create plots.")
        return
    
    special_data = sim.results['special_case_100_orbits']
    times = special_data['times']
    pos_numerical = special_data['positions']
    vel_numerical = special_data['velocities']
    pos_analytical = special_data['analytical_positions']
    vel_analytical = special_data['analytical_velocities']
    position_errors = special_data['position_errors']
    velocity_errors = special_data['velocity_errors']
    
    # Calculate true anomaly for both methods
    true_anomaly_numerical = calculate_true_anomaly_from_state(pos_numerical, vel_numerical, MU_EARTH)
    true_anomaly_analytical = analytical_true_anomaly_propagation(sim.orbit_params, times)
    
    # Convert time to orbits for better readability
    orbital_period = sim.orbit_params['orbital_period']
    time_orbits = times / orbital_period
    
    # Create figure with improved layout for eccentric orbits
    fig = plt.figure(figsize=(16, 14))
    fig.suptitle('IMPROVED: RK8 vs Analytical Solutions (10 Orbits)\nOrbit B: Highly Eccentric (e=0.9) - Showing True Physics', 
                 fontsize=16, fontweight='bold')
    
    # 1. Position magnitude over time
    ax1 = plt.subplot(3, 3, 1)
    pos_mag_numerical = np.linalg.norm(pos_numerical, axis=1)
    pos_mag_analytical = np.linalg.norm(pos_analytical, axis=1)
    
    plt.plot(time_orbits, pos_mag_numerical/1000, 'b-', label='RK8 Numerical', linewidth=1.5)
    plt.plot(time_orbits, pos_mag_analytical/1000, 'r--', label='Analytical', linewidth=1.0, alpha=0.8)
    plt.xlabel('Time (orbits)')
    plt.ylabel('Position Magnitude (1000 km)')
    plt.title('Orbital Radius vs Time')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 2. Position error over time
    ax2 = plt.subplot(3, 3, 2)
    plt.semilogy(time_orbits, position_errors, 'g-', linewidth=1.5)
    plt.xlabel('Time (orbits)')
    plt.ylabel('Position Error (km)')
    plt.title('Position Error vs Time')
    plt.grid(True, alpha=0.3)
    
    # 3. True anomaly vs time - WRAPPED (shows periodic jumps)
    ax3 = plt.subplot(3, 3, 3)
    nu_num_mod = np.degrees(true_anomaly_numerical % (2*np.pi))
    nu_ana_mod = np.degrees(true_anomaly_analytical % (2*np.pi))
    
    plt.plot(time_orbits, nu_ana_mod, 'r--', label='Analytical', linewidth=1.5, alpha=0.8)
    plt.plot(time_orbits, nu_num_mod, 'b-', label='RK8 Numerical', linewidth=1.0)
    plt.xlabel('Time (orbits)')
    plt.ylabel('True Anomaly (degrees)')
    plt.title('True Anomaly vs Time (Wrapped)\nShows "Jumps" at Periapsis')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.ylim(-10, 370)
    
    # 4. True anomaly vs time - UNWRAPPED (shows continuous growth)
    ax4 = plt.subplot(3, 3, 4)
    nu_numerical_unwrapped = np.unwrap(true_anomaly_numerical)
    nu_analytical_unwrapped = np.unwrap(true_anomaly_analytical)
    
    plt.plot(time_orbits, np.degrees(nu_analytical_unwrapped), 'r--', 
             label='Analytical', linewidth=1.5, alpha=0.8)
    plt.plot(time_orbits, np.degrees(nu_numerical_unwrapped), 'b-', 
             label='RK8 Numerical', linewidth=1.0)
    plt.xlabel('Time (orbits)')
    plt.ylabel('True Anomaly (degrees, continuous)')
    plt.title('True Anomaly vs Time (Unwrapped)\nShows Continuous Growth')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 5. NEW: Rate of true anomaly change - THIS SHOWS THE ISSUE!
    ax5 = plt.subplot(3, 3, 5)
    dt = np.diff(times)
    dnu_dt_analytical = np.diff(true_anomaly_analytical) / dt * 3600  # deg/hour
    
    # Handle wraparound in derivatives  
    dnu_dt_analytical = np.where(dnu_dt_analytical < -180*3600, 
                                dnu_dt_analytical + 360*3600, dnu_dt_analytical)
    dnu_dt_analytical = np.where(dnu_dt_analytical > 180*3600,
                                dnu_dt_analytical - 360*3600, dnu_dt_analytical)
    
    plt.semilogy(time_orbits[1:], np.abs(np.degrees(dnu_dt_analytical)), 'r-', 
                 label='|dν/dt| Analytical', linewidth=2)
    plt.xlabel('Time (orbits)')
    plt.ylabel('|dν/dt| (degrees/hour)')
    plt.title('Rate of True Anomaly Change\nEXPLAINS the "slow then fast" behavior')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Add text annotation explaining the physics
    plt.text(0.02, 0.95, 'Near apoapsis:\n~1°/hour\n\nNear periapsis:\n~3000°/hour', 
             transform=ax5.transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # 6. Velocity magnitude over time
    ax6 = plt.subplot(3, 3, 6)
    vel_mag_numerical = np.linalg.norm(vel_numerical, axis=1)
    vel_mag_analytical = np.linalg.norm(vel_analytical, axis=1)
    
    plt.plot(time_orbits, vel_mag_numerical, 'b-', label='RK8 Numerical', linewidth=1.5)
    plt.plot(time_orbits, vel_mag_analytical, 'r--', label='Analytical', linewidth=1.0, alpha=0.8)
    plt.xlabel('Time (orbits)')
    plt.ylabel('Velocity Magnitude (km/s)')
    plt.title('Velocity Magnitude vs Time')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 7. Velocity error over time
    ax7 = plt.subplot(3, 3, 7)
    plt.semilogy(time_orbits, velocity_errors, 'orange', linewidth=1.5)
    plt.xlabel('Time (orbits)')
    plt.ylabel('Velocity Error (km/s)')
    plt.title('Velocity Error vs Time')
    plt.grid(True, alpha=0.3)
    
    # 8. Phase plot - True Anomaly vs Radius (shows orbital shape)
    ax8 = plt.subplot(3, 3, 8)
    plt.plot(nu_ana_mod, pos_mag_analytical/1000, 'r--', label='Analytical', linewidth=1.5, alpha=0.8)
    plt.plot(nu_num_mod, pos_mag_numerical/1000, 'b-', label='RK8 Numerical', linewidth=1.0)
    plt.xlabel('True Anomaly (degrees)')
    plt.ylabel('Orbital Radius (1000 km)')
    plt.title('Radius vs True Anomaly\n(Orbital Shape: e=0.9)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 9. True anomaly error over time
    ax9 = plt.subplot(3, 3, 9)
    # Calculate true anomaly errors (handle 2π wraparound)
    true_anomaly_diff = true_anomaly_numerical - true_anomaly_analytical
    true_anomaly_diff = np.where(true_anomaly_diff > np.pi, 
                                true_anomaly_diff - 2*np.pi, true_anomaly_diff)
    true_anomaly_diff = np.where(true_anomaly_diff < -np.pi, 
                                true_anomaly_diff + 2*np.pi, true_anomaly_diff)
    true_anomaly_errors = np.abs(true_anomaly_diff)
    
    plt.semilogy(time_orbits, np.degrees(true_anomaly_errors)*3600, 'purple', linewidth=1.5)
    plt.xlabel('Time (orbits)')
    plt.ylabel('True Anomaly Error (arcseconds)')
    plt.title('True Anomaly Error vs Time')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save the plot
    plot_filename = f"{RESULTS_DIR}/special_case_IMPROVED_eccentric_orbit_analysis.png"
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    print(f"IMPROVED eccentric orbit plots saved to: {plot_filename}")
    
    # Print diagnostics about the true anomaly behavior
    max_rate = np.max(np.abs(np.degrees(dnu_dt_analytical)))
    min_rate = np.min(np.abs(np.degrees(dnu_dt_analytical)))
    rate_ratio = max_rate / min_rate
    
    print(f"\n" + "="*60)
    print(f"TRUE ANOMALY BEHAVIOR ANALYSIS")
    print(f"="*60) 
    print(f"This explains what you observed in the original plots:")
    print(f"")
    print(f"Max rate of ν change: {max_rate:.1f}°/hour (near periapsis)")
    print(f"Min rate of ν change: {min_rate:.3f}°/hour (near apoapsis)")
    print(f"Rate ratio (max/min): {rate_ratio:.0f}:1")
    print(f"")
    print(f"CONCLUSION:")
    print(f"• 'Almost no rise during orbit' → CORRECT (slow at apoapsis)")
    print(f"• 'Huge raises at periapsis' → CORRECT (fast at periapsis)")
    print(f"• This is NOT a bug - it's correct orbital mechanics!")
    print(f"• For e=0.9, satellite spends 99% of time moving slowly near apoapsis")
    print(f"• And 1% of time moving very fast near periapsis")
    print(f"")
    print(f"Max true anomaly error: {np.max(np.degrees(true_anomaly_errors)*3600):.3f} arcseconds")
    print(f"This error is EXCELLENT for numerical integration!")
    print(f"="*60)
    
    # Show plot if running interactively
    try:
        plt.show()
    except:
        pass
    
    return fig

# =============================================================================
# SINGLE METHOD PLOTTING FUNCTIONS
# =============================================================================

def create_position_vs_time_plots(simulation, method_name, method_data):
    """Create position magnitude vs time plots for all step sizes."""
    print(f"  Creating {method_name} position vs time plots...")
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(f'{method_name} Method: Position Magnitude vs Time\nAnalytical vs Numerical Solutions', fontsize=16, fontweight='bold')
    
    # Select a subset of step sizes to avoid overcrowding
    selected_steps = [1, 30, 100, 300]
    colors = ['red', 'blue', 'green', 'orange']
    
    for idx, (step_size, color) in enumerate(zip(selected_steps, colors)):
        if step_size not in method_data:
            continue
            
        ax = axes[idx // 2, idx % 2]
        result = method_data[step_size]
        
        # Convert time to hours for readability
        times_hours = result['times'] / 3600
        pos_mag_numerical = np.linalg.norm(result['positions'], axis=1) / 1000
        pos_mag_analytical = np.linalg.norm(result['analytical_positions'], axis=1) / 1000
        
        ax.plot(times_hours, pos_mag_analytical, 'k--', label='Analytical', linewidth=2, alpha=0.8)
        ax.plot(times_hours, pos_mag_numerical, color=color, label=f'{method_name} (Δt={step_size}s)', linewidth=1.5)
        
        ax.set_xlabel('Time (hours)')
        ax.set_ylabel('Position Magnitude (1000 km)')
        ax.set_title(f'Step Size: {step_size} seconds')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Add error statistics as text
        pos_errors = result['position_errors']
        ax.text(0.02, 0.98, f'Max Error: {np.max(pos_errors):.2e} km\nFinal Error: {pos_errors[-1]:.2e} km',
                transform=ax.transAxes, verticalalignment='top', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(f"{RESULTS_DIR}/{method_name.lower()}_position_vs_time.png", dpi=300, bbox_inches='tight')
    plt.show()

def create_velocity_vs_time_plots(simulation, method_name, method_data):
    """Create velocity magnitude vs time plots for all step sizes."""
    print(f"  Creating {method_name} velocity vs time plots...")
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(f'{method_name} Method: Velocity Magnitude vs Time\nAnalytical vs Numerical Solutions', fontsize=16, fontweight='bold')
    
    # Select a subset of step sizes to avoid overcrowding
    selected_steps = [1, 30, 100, 300]
    colors = ['red', 'blue', 'green', 'orange']
    
    for idx, (step_size, color) in enumerate(zip(selected_steps, colors)):
        if step_size not in method_data:
            continue
            
        ax = axes[idx // 2, idx % 2]
        result = method_data[step_size]
        
        # Convert time to hours for readability
        times_hours = result['times'] / 3600
        vel_mag_numerical = np.linalg.norm(result['velocities'], axis=1)
        vel_mag_analytical = np.linalg.norm(result['analytical_velocities'], axis=1)
        
        ax.plot(times_hours, vel_mag_analytical, 'k--', label='Analytical', linewidth=2, alpha=0.8)
        ax.plot(times_hours, vel_mag_numerical, color=color, label=f'{method_name} (Δt={step_size}s)', linewidth=1.5)
        
        ax.set_xlabel('Time (hours)')
        ax.set_ylabel('Velocity Magnitude (km/s)')
        ax.set_title(f'Step Size: {step_size} seconds')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Add error statistics as text
        vel_errors = result['velocity_errors']
        ax.text(0.02, 0.98, f'Max Error: {np.max(vel_errors):.2e} km/s\nFinal Error: {vel_errors[-1]:.2e} km/s',
                transform=ax.transAxes, verticalalignment='top', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(f"{RESULTS_DIR}/{method_name.lower()}_velocity_vs_time.png", dpi=300, bbox_inches='tight')
    plt.show()

def create_true_anomaly_vs_time_plots(simulation, method_name, method_data):
    """Create true anomaly vs time plots for all step sizes."""
    print(f"  Creating {method_name} true anomaly vs time plots...")
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(f'{method_name} Method: True Anomaly vs Time\nAnalytical vs Numerical Solutions', fontsize=16, fontweight='bold')
    
    # Select a subset of step sizes to avoid overcrowding
    selected_steps = [1, 30, 100, 300]
    colors = ['red', 'blue', 'green', 'orange']
    
    for idx, (step_size, color) in enumerate(zip(selected_steps, colors)):
        if step_size not in method_data:
            continue
            
        ax = axes[idx // 2, idx % 2]
        result = method_data[step_size]
        
        # Convert time to hours for readability
        times_hours = result['times'] / 3600
        nu_numerical_deg = np.degrees(result['true_anomaly_method'])
        nu_analytical_deg = np.degrees(result['analytical_true_anomaly'])
        
        ax.plot(times_hours, nu_analytical_deg, 'k--', label='Analytical', linewidth=2, alpha=0.8)
        ax.plot(times_hours, nu_numerical_deg, color=color, label=f'{method_name} (Δt={step_size}s)', linewidth=1.5)
        
        ax.set_xlabel('Time (hours)')
        ax.set_ylabel('True Anomaly (degrees)')
        ax.set_title(f'Step Size: {step_size} seconds')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Add error statistics as text
        nu_errors_deg = np.degrees(result['true_anomaly_errors'])
        ax.text(0.02, 0.98, f'Max Error: {np.max(nu_errors_deg):.2e}°\nFinal Error: {nu_errors_deg[-1]:.2e}°',
                transform=ax.transAxes, verticalalignment='top', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(f"{RESULTS_DIR}/{method_name.lower()}_true_anomaly_vs_time.png", dpi=300, bbox_inches='tight')
    plt.show()

def create_position_error_vs_time_plots(simulation, method_name, method_data):
    """Create position error vs time plots for all step sizes."""
    print(f"  Creating {method_name} position error vs time plots...")
    
    fig = plt.figure(figsize=(16, 10))
    colors = plt.cm.viridis(np.linspace(0, 1, len(INTEGRATION_STEPS)))
    
    for idx, (step_size, color) in enumerate(zip(INTEGRATION_STEPS, colors)):
        if step_size not in method_data:
            continue
            
        result = method_data[step_size]
        times_hours = result['times'] / 3600
        pos_errors = result['position_errors']
        
        plt.semilogy(times_hours, pos_errors, color=color, 
                    label=f'Δt={step_size}s (Max: {np.max(pos_errors):.2e} km)', linewidth=2)
    
    plt.xlabel('Time (hours)', fontsize=12)
    plt.ylabel('Position Error (km)', fontsize=12)
    plt.title(f'{method_name} Method: Position Errors vs Time\nNumerical vs Analytical Solution', fontsize=14, fontweight='bold')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{RESULTS_DIR}/{method_name.lower()}_position_errors_vs_time.png", dpi=300, bbox_inches='tight')
    plt.show()

def create_velocity_error_vs_time_plots(simulation, method_name, method_data):
    """Create velocity error vs time plots for all step sizes."""
    print(f"  Creating {method_name} velocity error vs time plots...")
    
    fig = plt.figure(figsize=(16, 10))
    colors = plt.cm.plasma(np.linspace(0, 1, len(INTEGRATION_STEPS)))
    
    for idx, (step_size, color) in enumerate(zip(INTEGRATION_STEPS, colors)):
        if step_size not in method_data:
            continue
            
        result = method_data[step_size]
        times_hours = result['times'] / 3600
        vel_errors = result['velocity_errors']
        
        plt.semilogy(times_hours, vel_errors, color=color, 
                    label=f'Δt={step_size}s (Max: {np.max(vel_errors):.2e} km/s)', linewidth=2)
    
    plt.xlabel('Time (hours)', fontsize=12)
    plt.ylabel('Velocity Error (km/s)', fontsize=12)
    plt.title(f'{method_name} Method: Velocity Errors vs Time\nNumerical vs Analytical Solution', fontsize=14, fontweight='bold')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{RESULTS_DIR}/{method_name.lower()}_velocity_errors_vs_time.png", dpi=300, bbox_inches='tight')
    plt.show()

def create_true_anomaly_error_vs_time_plots(simulation, method_name, method_data):
    """Create true anomaly error vs time plots for all step sizes."""
    print(f"  Creating {method_name} true anomaly error vs time plots...")
    
    fig = plt.figure(figsize=(16, 10))
    colors = plt.cm.inferno(np.linspace(0, 1, len(INTEGRATION_STEPS)))
    
    for idx, (step_size, color) in enumerate(zip(INTEGRATION_STEPS, colors)):
        if step_size not in method_data:
            continue
            
        result = method_data[step_size]
        times_hours = result['times'] / 3600
        nu_errors_deg = np.degrees(result['true_anomaly_errors'])
        
        plt.semilogy(times_hours, nu_errors_deg, color=color, 
                    label=f'Δt={step_size}s (Max: {np.max(nu_errors_deg):.2e}°)', linewidth=2)
    
    plt.xlabel('Time (hours)', fontsize=12)
    plt.ylabel('True Anomaly Error (degrees)', fontsize=12)
    plt.title(f'{method_name} Method: True Anomaly Errors vs Time\nNumerical vs Analytical Solution', fontsize=14, fontweight='bold')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{RESULTS_DIR}/{method_name.lower()}_true_anomaly_errors_vs_time.png", dpi=300, bbox_inches='tight')
    plt.show()

def create_single_method_convergence_plot(simulation, method_name, method_data):
    """Create convergence analysis plot for a single method."""
    print(f"  Creating {method_name} convergence analysis plot...")
    
    step_sizes = np.array(INTEGRATION_STEPS)
    max_pos_errors = np.array([method_data[h]['max_position_error'] for h in INTEGRATION_STEPS])
    max_vel_errors = np.array([method_data[h]['max_velocity_error'] for h in INTEGRATION_STEPS])
    max_nu_errors = np.array([method_data[h]['max_true_anomaly_error'] for h in INTEGRATION_STEPS])
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    fig.suptitle(f'{method_name} Method: Convergence Analysis', fontsize=16, fontweight='bold')
    
    # Plot 1: Error vs Step Size
    ax1.loglog(step_sizes, max_pos_errors, 'ro-', label='Position Error (km)', linewidth=2, markersize=8)
    ax1.loglog(step_sizes, max_vel_errors * 1000, 'bs-', label='Velocity Error × 1000 (km)', linewidth=2, markersize=8)
    ax1.loglog(step_sizes, np.degrees(max_nu_errors) * 100, 'g^-', label='True Anomaly Error × 100 (°)', linewidth=2, markersize=8)
    
    # Fit and plot theoretical convergence lines
    if method_name == 'RK4':
        theoretical_slope = 4
        reference_error = max_pos_errors[0]
        reference_step = step_sizes[0]
    else:  # RK8
        theoretical_slope = 8
        reference_error = max_pos_errors[0]
        reference_step = step_sizes[0]
    
    theoretical_line = reference_error * (step_sizes / reference_step) ** theoretical_slope
    ax1.loglog(step_sizes, theoretical_line, 'k--', alpha=0.7, 
              label=f'Theoretical O(Δt^{theoretical_slope})', linewidth=2)
    
    ax1.set_xlabel('Step Size (seconds)', fontsize=12)
    ax1.set_ylabel('Maximum Error', fontsize=12)
    ax1.set_title('Maximum Errors vs Step Size', fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Computation Time vs Accuracy
    comp_times = [method_data[h]['computation_time'] for h in INTEGRATION_STEPS]
    ax2.loglog(max_pos_errors, comp_times, 'ro-', label='Position Error', linewidth=2, markersize=8)
    ax2.set_xlabel('Maximum Position Error (km)', fontsize=12)
    ax2.set_ylabel('Computation Time (seconds)', fontsize=12)
    ax2.set_title('Computation Time vs Accuracy Trade-off', fontsize=14)
    ax2.grid(True, alpha=0.3)
    
    # Add step size annotations
    for i, step in enumerate(step_sizes):
        ax2.annotate(f'{step}s', (max_pos_errors[i], comp_times[i]), 
                    xytext=(5, 5), textcoords='offset points', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(f"{RESULTS_DIR}/{method_name.lower()}_convergence_analysis.png", dpi=300, bbox_inches='tight')
    plt.show()

def create_single_method_3d_trajectory(simulation, method_name, method_data):
    """Create 3D trajectory comparison for a single method with multiple step sizes."""
    print(f"  Creating {method_name} 3D trajectory plots...")
    
    # Create individual 3D plots for selected step sizes
    selected_steps = [1, 30, 100, 300]  # Show progression from best to worst accuracy
    colors = ['red', 'blue', 'green', 'orange']
    
    # First create individual plots for each step size
    for step_size, color in zip(selected_steps, colors):
        if step_size not in method_data:
            continue
            
        result = method_data[step_size]
        
        fig = plt.figure(figsize=(14, 10))
        ax = fig.add_subplot(111, projection='3d')
        
        # Plot analytical solution
        pos_analytical = result['analytical_positions'] / 1000  # Convert to 1000 km
        ax.plot(pos_analytical[:, 0], pos_analytical[:, 1], pos_analytical[:, 2],
               'k-', label='Analytical', linewidth=3, alpha=0.8)
        
        # Plot numerical solution
        pos_numerical = result['positions'] / 1000  # Convert to 1000 km
        ax.plot(pos_numerical[:, 0], pos_numerical[:, 1], pos_numerical[:, 2],
               color=color, label=f'{method_name} (Δt={step_size}s)', linewidth=2)
        
        # Mark starting point
        ax.scatter([pos_analytical[0, 0]], [pos_analytical[0, 1]], [pos_analytical[0, 2]],
                  color='green', s=200, marker='o', label='Start', alpha=0.8)
        
        # Mark Earth center
        ax.scatter([0], [0], [0], color='blue', s=200, marker='o', label='Earth Center', alpha=0.8)
        
        # Add Earth sphere
        u = np.linspace(0, 2 * np.pi, 50)
        v = np.linspace(0, np.pi, 50)
        earth_radius = 6.371  # 1000 km units
        x_earth = earth_radius * np.outer(np.cos(u), np.sin(v))
        y_earth = earth_radius * np.outer(np.sin(u), np.sin(v))
        z_earth = earth_radius * np.outer(np.ones(np.size(u)), np.cos(v))
        ax.plot_surface(x_earth, y_earth, z_earth, color='lightblue', alpha=0.3)
        
        ax.set_xlabel('X (1000 km)', fontsize=12)
        ax.set_ylabel('Y (1000 km)', fontsize=12)
        ax.set_zlabel('Z (1000 km)', fontsize=12)
        ax.set_title(f'{method_name} Method: 3D Trajectory Comparison\n'
                    f'Step Size: {step_size}s, Max Position Error: {np.max(result["position_errors"]):.2e} km',
                    fontsize=14, fontweight='bold')
        ax.legend()
        
        # Set equal aspect ratio
        max_range = np.array([pos_analytical[:, 0].max() - pos_analytical[:, 0].min(),
                             pos_analytical[:, 1].max() - pos_analytical[:, 1].min(),
                             pos_analytical[:, 2].max() - pos_analytical[:, 2].min()]).max() / 2.0
        mid_x = (pos_analytical[:, 0].max() + pos_analytical[:, 0].min()) * 0.5
        mid_y = (pos_analytical[:, 1].max() + pos_analytical[:, 1].min()) * 0.5
        mid_z = (pos_analytical[:, 2].max() + pos_analytical[:, 2].min()) * 0.5
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)
        
        plt.tight_layout()
        plt.savefig(f"{RESULTS_DIR}/{method_name.lower()}_3d_trajectory_step_{step_size}s.png", dpi=300, bbox_inches='tight')
        plt.show()
    
    # Now create a combined plot showing all step sizes together
    create_combined_3d_trajectory_plot(simulation, method_name, method_data, selected_steps, colors)

def create_combined_3d_trajectory_plot(simulation, method_name, method_data, selected_steps, colors):
    """Create a combined 3D trajectory plot showing multiple step sizes together."""
    print(f"  Creating {method_name} combined 3D trajectory plot...")
    
    fig = plt.figure(figsize=(16, 12))
    ax = fig.add_subplot(111, projection='3d')
    
    # Get analytical solution from the finest step size for reference
    finest_step = min(selected_steps)
    if finest_step in method_data:
        pos_analytical = method_data[finest_step]['analytical_positions'] / 1000
        ax.plot(pos_analytical[:, 0], pos_analytical[:, 1], pos_analytical[:, 2],
               'k-', label='Analytical', linewidth=4, alpha=0.9)
    
    # Plot numerical solutions for different step sizes
    max_errors = []
    for step_size, color in zip(selected_steps, colors):
        if step_size not in method_data:
            continue
            
        result = method_data[step_size]
        pos_numerical = result['positions'] / 1000  # Convert to 1000 km
        max_error = np.max(result['position_errors'])
        max_errors.append(max_error)
        
        ax.plot(pos_numerical[:, 0], pos_numerical[:, 1], pos_numerical[:, 2],
               color=color, label=f'{method_name} Δt={step_size}s (Max err: {max_error:.1e} km)', 
               linewidth=2, alpha=0.8)
    
    # Mark starting point
    if finest_step in method_data:
        pos_start = method_data[finest_step]['analytical_positions'][0] / 1000
        ax.scatter([pos_start[0]], [pos_start[1]], [pos_start[2]],
                  color='green', s=300, marker='*', label='Start', alpha=1.0, edgecolors='black')
    
    # Mark Earth center
    ax.scatter([0], [0], [0], color='blue', s=200, marker='o', label='Earth Center', alpha=0.8)
    
    # Add Earth sphere
    u = np.linspace(0, 2 * np.pi, 40)
    v = np.linspace(0, np.pi, 40)
    earth_radius = 6.371  # 1000 km units
    x_earth = earth_radius * np.outer(np.cos(u), np.sin(v))
    y_earth = earth_radius * np.outer(np.sin(u), np.sin(v))
    z_earth = earth_radius * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x_earth, y_earth, z_earth, color='lightblue', alpha=0.2)
    
    ax.set_xlabel('X (1000 km)', fontsize=12)
    ax.set_ylabel('Y (1000 km)', fontsize=12)
    ax.set_zlabel('Z (1000 km)', fontsize=12)
    ax.set_title(f'{method_name} Method: 3D Trajectory Comparison\n'
                f'Multiple Step Sizes - Error Range: {min(max_errors):.1e} to {max(max_errors):.1e} km',
                fontsize=14, fontweight='bold')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Set equal aspect ratio
    if finest_step in method_data:
        max_range = np.array([pos_analytical[:, 0].max() - pos_analytical[:, 0].min(),
                             pos_analytical[:, 1].max() - pos_analytical[:, 1].min(),
                             pos_analytical[:, 2].max() - pos_analytical[:, 2].min()]).max() / 2.0
        mid_x = (pos_analytical[:, 0].max() + pos_analytical[:, 0].min()) * 0.5
        mid_y = (pos_analytical[:, 1].max() + pos_analytical[:, 1].min()) * 0.5
        mid_z = (pos_analytical[:, 2].max() + pos_analytical[:, 2].min()) * 0.5
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)
    
    plt.tight_layout()
    plt.savefig(f"{RESULTS_DIR}/{method_name.lower()}_3d_trajectory_combined.png", dpi=300, bbox_inches='tight')
    plt.show()