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

from core.Constants import INTEGRATION_STEPS, RESULTS_DIR

# =============================================================================
# MAIN PLOTTING FUNCTION
# =============================================================================

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