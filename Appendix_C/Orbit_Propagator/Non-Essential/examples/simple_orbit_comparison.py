"""
Simple Orbit Comparison: Analytical vs 100s RK4 Integration
===========================================================

This script creates a clean plot comparing the analytical orbit solution
with the RK4 numerical integration using 100s time steps.

Author: Generated for SMAD Moon Communications Constellation Study
Date: October 2025
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from pathlib import Path

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import our modules
from core.Constants import (
    get_orbit_b_parameters, 
    RESULTS_DIR,
    print_orbit_summary
)
from core.Kepler import analytical_propagation, position_velocity_from_elements
from core.RK4 import propagate_orbit_rk4

def plot_analytical_vs_100s():
    """
    Create a simple comparison plot of analytical vs 100s RK4 orbit.
    """
    print("Analytical vs 100s Time Step Orbit Comparison")
    print("=" * 50)
    
    # Get orbital parameters
    orbit_params = get_orbit_b_parameters()
    print_orbit_summary()
    
    # Simulation parameters
    analytical_step_size = 50.0   # Fine step for analytical (smooth curve)
    rk4_step_size = 300.0         # Larger step for RK4 (to show differences)
    num_orbits = 2                # Show 2 orbits to see error accumulation
    period = orbit_params['orbital_period']
    t_end = num_orbits * period
    
    print(f"\nSimulation Parameters:")
    print(f"Analytical step size: {analytical_step_size} s")
    print(f"RK4 step size: {rk4_step_size} s")
    print(f"Number of orbits: {num_orbits}")
    print(f"Simulation time: {t_end/3600:.2f} hours")
    
    # Get initial conditions
    a = orbit_params['semimajor_axis']
    e = orbit_params['eccentricity']
    i = orbit_params['inclination']
    raan = orbit_params['raan']
    arg_p = orbit_params['arg_periapsis']
    nu_0 = orbit_params['true_anomaly_0']
    
    r0, v0 = position_velocity_from_elements(a, e, i, raan, arg_p, nu_0)
    initial_state = np.concatenate([r0, v0])
    
    print(f"\nRunning simulations...")
    
    # Generate analytical solution with fine step size (smooth curve)
    time_analytical = np.arange(0, t_end + analytical_step_size, analytical_step_size)
    pos_analytical, vel_analytical = analytical_propagation(orbit_params, time_analytical)
    
    # Run RK4 integration with larger step size (to show numerical errors)
    times_rk4, pos_rk4, vel_rk4 = propagate_orbit_rk4(
        initial_state,
        (0, t_end),
        rk4_step_size
    )
    
    print(f"Analytical points: {len(time_analytical)}")
    print(f"RK4 points: {len(times_rk4)}")
    
    # Create the comparison plot
    create_orbit_comparison_plot(
        pos_analytical, pos_rk4, 
        time_analytical, times_rk4,
        orbit_params, rk4_step_size, num_orbits
    )
    
    # Create the perigee zoom plot
    create_perigee_zoom_plot(
        pos_analytical, pos_rk4, 
        time_analytical, times_rk4,
        orbit_params, rk4_step_size, num_orbits
    )

def create_orbit_comparison_plot(pos_analytical, pos_rk4, times_analytical, times_rk4,
                                orbit_params, step_size, num_orbits):
    """
    Create a clean 2D comparison plot of the two orbits.
    """
    print(f"Creating orbit comparison plot...")
    
    # Create figure with single 2D plot
    fig = plt.figure(figsize=(12, 10))
    
    # Define colors
    analytical_color = "#000000"  # Black
    rk4_orbit1_color = '#DC143C'  # Crimson (first orbit)
    rk4_orbit2_color = '#FF8C00'  # Dark Orange (second orbit)
    
    # 2D XY Plane View (clearer for highly eccentric orbit)
    ax = fig.add_subplot(1, 1, 1)
    
    # Calculate period points to separate orbits
    period = orbit_params['orbital_period']
    
    # Split data by orbits for different coloring
    # Find indices for first and second orbit
    orbit1_mask = times_rk4 <= period
    orbit2_mask = times_rk4 > period
    
    # Plot analytical orbit (all orbits in black)
    ax.plot(pos_analytical[:, 0]/1000, pos_analytical[:, 1]/1000,
             color=analytical_color, linewidth=1, label='Analytical Solution', alpha=0.8)
    
    # Plot RK4 first orbit
    if np.any(orbit1_mask):
        ax.plot(pos_rk4[orbit1_mask, 0]/1000, pos_rk4[orbit1_mask, 1]/1000,
                 color=rk4_orbit1_color, linewidth=1, linestyle='-', 
                 label=f'RK4 Orbit 1 (Δt={step_size}s)', alpha=0.9)
        
        # Add markers for first orbit
        rk4_points_1 = pos_rk4[orbit1_mask]
        every_nth = max(1, len(rk4_points_1) // 15)
        ax.scatter(rk4_points_1[::every_nth, 0]/1000, rk4_points_1[::every_nth, 1]/1000,
                   color=rk4_orbit1_color, s=8, alpha=0.6, zorder=5)
    
    # Plot RK4 second orbit  
    if np.any(orbit2_mask):
        ax.plot(pos_rk4[orbit2_mask, 0]/1000, pos_rk4[orbit2_mask, 1]/1000,
                 color=rk4_orbit2_color, linewidth=1, linestyle='-', 
                 label=f'RK4 Orbit 2 (Δt={step_size}s)', alpha=0.9)
        
        # Add markers for second orbit
        rk4_points_2 = pos_rk4[orbit2_mask]
        every_nth = max(1, len(rk4_points_2) // 15)
        ax.scatter(rk4_points_2[::every_nth, 0]/1000, rk4_points_2[::every_nth, 1]/1000,
                   color=rk4_orbit2_color, s=8, alpha=0.6, zorder=5)
    
    # Mark key points
    ax.scatter([orbit_params['r_periapsis']/1000], [0], 
               color='red', s=150, marker='*', label='Periapsis', zorder=10)
    ax.scatter([-orbit_params['r_apoapsis']/1000], [0], 
               color='blue', s=150, marker='*', label='Apoapsis', zorder=10)
    ax.scatter([0], [0], color='black', s=100, marker='o', 
               label='Earth Center', zorder=10)
    
    # Add Earth radius circle for reference
    earth_radius = 6371  # km
    circle = plt.Circle((0, 0), earth_radius/1000, color='lightblue', 
                       fill=False, linestyle=':', alpha=0.7, label='Earth Surface')
    ax.add_patch(circle)
    
    ax.set_xlabel('X (1000 km)', fontsize=14)
    ax.set_ylabel('Y (1000 km)', fontsize=14)
    ax.set_title(f'Orbit Comparison: Analytical vs RK4\\n{num_orbits} Orbits | RK4 Step Size: {step_size}s | e = {orbit_params["eccentricity"]:.3f}', 
                  fontsize=16, fontweight='bold')
    ax.legend(loc='upper right', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    
    # Add orbital parameters as text
    param_text = f'''Orbital Parameters:
Semi-major axis: {orbit_params['semimajor_axis']/1000:.1f} × 1000 km
Eccentricity: {orbit_params['eccentricity']:.3f}
Periapsis: {orbit_params['r_periapsis']/1000:.1f} × 1000 km
Apoapsis: {orbit_params['r_apoapsis']/1000:.1f} × 1000 km
Period: {orbit_params['orbital_period']/3600:.2f} hours'''
    
    ax.text(0.02, 0.98, param_text, transform=ax.transAxes, 
             verticalalignment='top', fontsize=11, 
             bbox=dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.8))
    
    plt.tight_layout()
    
    # Save the plot
    results_dir = Path(RESULTS_DIR)
    results_dir.mkdir(exist_ok=True)
    
    filename = f'analytical_vs_100s_orbit_comparison.png'
    filepath = results_dir / filename
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    
    print(f"Plot saved: {filepath}")
    plt.show()

def create_perigee_zoom_plot(pos_analytical, pos_rk4, times_analytical, times_rk4,
                            orbit_params, step_size, num_orbits):
    """
    Create a zoomed-in plot of the perigee region to show detailed behavior.
    """
    print(f"Creating perigee zoom plot...")
    
    # Create figure
    fig = plt.figure(figsize=(12, 10))
    
    # Define colors
    analytical_color = '#2E8B57'  # Sea Green
    rk4_color = '#DC143C'        # Crimson
    
    # 2D XY Plane View
    ax = fig.add_subplot(1, 1, 1)
    
    # Plot analytical orbit
    ax.plot(pos_analytical[:, 0]/1000, pos_analytical[:, 1]/1000,
             color=analytical_color, linewidth=1, label='Analytical Solution', alpha=0.8)
    
    # Plot RK4 orbit
    ax.plot(pos_rk4[:, 0]/1000, pos_rk4[:, 1]/1000,
             color=rk4_color, linewidth=1, linestyle='-', 
             label=f'RK4 (Δt={step_size}s)', alpha=0.9)
    
    # Add markers for RK4 points
    every_nth = max(1, len(pos_rk4) // 30)  # Show more markers in zoom
    ax.scatter(pos_rk4[::every_nth, 0]/1000, pos_rk4[::every_nth, 1]/1000,
               color=rk4_color, s=12, alpha=0.6, zorder=5)
    
    # Mark perigee point
    perigee_x = orbit_params['r_periapsis']/1000
    ax.scatter([perigee_x], [0], 
               color='red', s=200, marker='*', label='Perigee', zorder=10)
    
    # Mark Earth center
    ax.scatter([0], [0], color='black', s=120, marker='o', 
               label='Earth Center', zorder=10)
    
    # Add Earth radius circle for reference
    earth_radius = 6371  # km
    circle = plt.Circle((0, 0), earth_radius/1000, color='lightblue', 
                       fill=False, linestyle=':', alpha=0.7, label='Earth Surface')
    ax.add_patch(circle)
    
    # Set zoom limits around perigee
    perigee_radius = orbit_params['r_periapsis']/1000
    zoom_factor = 1.01  # Show region 1.01x the perigee radius (much more zoomed)
    zoom_range = perigee_radius * zoom_factor
    
    ax.set_xlim(-zoom_range * 0.3, zoom_range * 1.2)
    ax.set_ylim(-zoom_range * 0.8, zoom_range * 0.8)
    
    ax.set_xlabel('X (1000 km)', fontsize=14)
    ax.set_ylabel('Y (1000 km)', fontsize=14)
    ax.set_title(f'Perigee Region Zoom: Analytical vs RK4\\nStep Size: {step_size}s | Perigee: {perigee_radius:.1f} × 1000 km', 
                  fontsize=16, fontweight='bold')
    ax.legend(loc='upper right', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    
    # Add zoom info as text
    zoom_text = f'''Zoom Details:
Perigee altitude: {(orbit_params['r_periapsis'] - earth_radius*1000)/1000:.1f} × 1000 km
Earth radius: {earth_radius/1000:.1f} × 1000 km
Zoom factor: {zoom_factor:.1f}x perigee radius
Eccentricity: {orbit_params['eccentricity']:.3f}'''
    
    ax.text(0.02, 0.98, zoom_text, transform=ax.transAxes, 
             verticalalignment='top', fontsize=11, 
             bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow', alpha=0.8))
    
    plt.tight_layout()
    
    # Save the plot
    results_dir = Path(RESULTS_DIR)
    results_dir.mkdir(exist_ok=True)
    
    filename = f'perigee_zoom_comparison.png'
    filepath = results_dir / filename
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    
    print(f"Perigee zoom plot saved: {filepath}")
    plt.show()

if __name__ == "__main__":
    plot_analytical_vs_100s()