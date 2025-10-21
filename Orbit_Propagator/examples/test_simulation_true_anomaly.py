#!/usr/bin/env python3
"""
Test the true anomaly behavior in actual orbital simulation
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.Constants import get_orbit_b_parameters, MU_EARTH
from core.Kepler import analytical_propagation, analytical_true_anomaly_propagation
from core.RK8 import propagate_orbit_rk8
from core.RK4 import calculate_true_anomaly_from_state
from core.Kepler import position_velocity_from_elements

def test_true_anomaly_in_simulation():
    """Test true anomaly calculation in a short simulation."""
    print("Testing True Anomaly in Orbital Simulation")
    print("=" * 50)
    
    # Get Orbit B parameters  
    orbit_params = get_orbit_b_parameters()
    print(f"Orbit parameters: a={orbit_params['semimajor_axis']} km, e={orbit_params['eccentricity']}")
    
    # Set up initial conditions
    a = orbit_params['semimajor_axis']
    e = orbit_params['eccentricity'] 
    i = orbit_params['inclination']
    raan = orbit_params['raan']
    arg_p = orbit_params['arg_periapsis']
    nu_0 = orbit_params['true_anomaly_0']
    
    r0, v0 = position_velocity_from_elements(a, e, i, raan, arg_p, nu_0)
    initial_state = np.concatenate([r0, v0])
    
    # Short simulation: 1/4 orbit around periapsis passage
    period = orbit_params['orbital_period']
    t_end = period / 4  # Quarter orbit
    step_size = period / 1000  # Small steps for resolution
    
    print(f"Simulation time: {t_end/3600:.2f} hours ({t_end/period:.3f} orbits)")
    print(f"Step size: {step_size:.1f} seconds")
    
    # Run RK8 integration
    times_rk8, pos_rk8, vel_rk8 = propagate_orbit_rk8(initial_state, (0, t_end), step_size)
    
    # Get analytical solution
    pos_analytical, vel_analytical = analytical_propagation(orbit_params, times_rk8)
    
    # Calculate true anomalies
    nu_rk8 = calculate_true_anomaly_from_state(pos_rk8, vel_rk8, MU_EARTH)
    nu_analytical = analytical_true_anomaly_propagation(orbit_params, times_rk8)
    
    # Convert times to orbital periods
    time_orbits = times_rk8 / period
    
    print(f"\nResults:")
    print(f"Time points: {len(times_rk8)}")
    print(f"True anomaly range (analytical): {np.degrees(np.min(nu_analytical)):.1f}° to {np.degrees(np.max(nu_analytical)):.1f}°")
    print(f"True anomaly range (numerical): {np.degrees(np.min(nu_rk8)):.1f}° to {np.degrees(np.max(nu_rk8)):.1f}°")
    
    # Check for discontinuities
    nu_analytical_diff = np.diff(nu_analytical)
    nu_rk8_diff = np.diff(nu_rk8)
    
    # Find large jumps (potential discontinuities)
    large_jumps_analytical = np.where(np.abs(nu_analytical_diff) > np.pi/2)[0]
    large_jumps_rk8 = np.where(np.abs(nu_rk8_diff) > np.pi/2)[0]
    
    print(f"\nLarge jumps in analytical ν: {len(large_jumps_analytical)}")
    if len(large_jumps_analytical) > 0:
        for idx in large_jumps_analytical[:5]:  # Show first 5
            print(f"  At t={time_orbits[idx]:.4f} orbits: Δν = {np.degrees(nu_analytical_diff[idx]):.1f}°")
    
    print(f"Large jumps in numerical ν: {len(large_jumps_rk8)}")  
    if len(large_jumps_rk8) > 0:
        for idx in large_jumps_rk8[:5]:  # Show first 5
            print(f"  At t={time_orbits[idx]:.4f} orbits: Δν = {np.degrees(nu_rk8_diff[idx]):.1f}°")
    
    # Create plots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
    
    # Plot 1: True anomaly vs time (raw)
    ax1.plot(time_orbits, np.degrees(nu_analytical), 'r-', label='Analytical', linewidth=2)
    ax1.plot(time_orbits, np.degrees(nu_rk8), 'b--', label='RK8 Numerical', alpha=0.7)
    ax1.set_xlabel('Time (orbits)')
    ax1.set_ylabel('True Anomaly (degrees)')
    ax1.set_title('True Anomaly vs Time (Raw)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: True anomaly vs time (unwrapped)
    nu_analytical_unwrapped = np.unwrap(nu_analytical)
    nu_rk8_unwrapped = np.unwrap(nu_rk8)
    ax2.plot(time_orbits, np.degrees(nu_analytical_unwrapped), 'r-', label='Analytical', linewidth=2)
    ax2.plot(time_orbits, np.degrees(nu_rk8_unwrapped), 'b--', label='RK8 Numerical', alpha=0.7)
    ax2.set_xlabel('Time (orbits)')
    ax2.set_ylabel('True Anomaly (degrees, unwrapped)')
    ax2.set_title('True Anomaly vs Time (Unwrapped)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: True anomaly differences
    ax3.plot(time_orbits[1:], np.degrees(nu_analytical_diff), 'r-', label='Analytical Δν', linewidth=2)
    ax3.plot(time_orbits[1:], np.degrees(nu_rk8_diff), 'b--', label='RK8 Numerical Δν', alpha=0.7)
    ax3.set_xlabel('Time (orbits)')
    ax3.set_ylabel('Δν per step (degrees)')
    ax3.set_title('True Anomaly Rate of Change')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Errors
    nu_error = nu_rk8 - nu_analytical
    # Handle wraparound in error calculation
    nu_error = np.where(nu_error > np.pi, nu_error - 2*np.pi, nu_error)
    nu_error = np.where(nu_error < -np.pi, nu_error + 2*np.pi, nu_error)
    
    ax4.plot(time_orbits, np.degrees(nu_error), 'g-', linewidth=2)
    ax4.set_xlabel('Time (orbits)')
    ax4.set_ylabel('True Anomaly Error (degrees)')
    ax4.set_title('RK8 vs Analytical True Anomaly Error')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('true_anomaly_analysis.png', dpi=300, bbox_inches='tight')
    print(f"\nPlots saved as 'true_anomaly_analysis.png'")
    
    # Show some statistics
    print(f"\nStatistics:")
    print(f"Max true anomaly error: {np.degrees(np.max(np.abs(nu_error))):.6f}°")
    print(f"RMS true anomaly error: {np.degrees(np.sqrt(np.mean(nu_error**2))):.6f}°")
    
    return time_orbits, nu_analytical, nu_rk8, nu_error

if __name__ == "__main__":
    test_true_anomaly_in_simulation()