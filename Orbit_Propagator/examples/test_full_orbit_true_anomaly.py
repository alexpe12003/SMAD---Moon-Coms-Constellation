#!/usr/bin/env python3
"""
Test true anomaly behavior across complete orbits with periapsis passage
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

def test_full_orbit_true_anomaly():
    """Test true anomaly across multiple complete orbits."""
    print("Testing True Anomaly Across Complete Orbits")
    print("=" * 50)
    
    # Get Orbit B parameters
    orbit_params = get_orbit_b_parameters()
    print(f"Orbit: a={orbit_params['semimajor_axis']} km, e={orbit_params['eccentricity']}")
    
    # Set up initial conditions
    a = orbit_params['semimajor_axis']
    e = orbit_params['eccentricity']
    i = orbit_params['inclination'] 
    raan = orbit_params['raan']
    arg_p = orbit_params['arg_periapsis']
    nu_0 = orbit_params['true_anomaly_0']
    
    r0, v0 = position_velocity_from_elements(a, e, i, raan, arg_p, nu_0)
    initial_state = np.concatenate([r0, v0])
    
    # Simulate 2 complete orbits
    period = orbit_params['orbital_period']
    t_end = 2 * period  # 2 complete orbits
    step_size = period / 500  # 500 steps per orbit
    
    print(f"Simulation: {t_end/period:.1f} orbits ({t_end/3600:.1f} hours)")
    print(f"Step size: {step_size:.1f} s ({500} steps per orbit)")
    print(f"Starting true anomaly: {np.degrees(nu_0):.1f}°")
    
    # Run RK8 integration
    times_rk8, pos_rk8, vel_rk8 = propagate_orbit_rk8(initial_state, (0, t_end), step_size)
    
    # Get analytical solution
    pos_analytical, vel_analytical = analytical_propagation(orbit_params, times_rk8)
    
    # Calculate true anomalies
    nu_rk8 = calculate_true_anomaly_from_state(pos_rk8, vel_rk8, MU_EARTH)
    nu_analytical = analytical_true_anomaly_propagation(orbit_params, times_rk8)
    
    # Convert times to orbital periods  
    time_orbits = times_rk8 / period
    
    print(f"\nSimulation completed:")
    print(f"Time points: {len(times_rk8)}")
    print(f"Final true anomaly (analytical): {np.degrees(nu_analytical[-1]):.1f}°")
    print(f"Final true anomaly (numerical): {np.degrees(nu_rk8[-1]):.1f}°")
    
    # Analyze true anomaly progression
    print(f"\nTrue Anomaly Analysis:")
    
    # Find periapsis passages (where true anomaly wraps from ~360° to ~0°)
    nu_analytical_deg = np.degrees(nu_analytical)
    nu_rk8_deg = np.degrees(nu_rk8)
    
    # Look for wraparound points (large negative jumps)
    nu_analytical_diff = np.diff(nu_analytical_deg)
    nu_rk8_diff = np.diff(nu_rk8_deg)
    
    wraparound_analytical = np.where(nu_analytical_diff < -180)[0]
    wraparound_rk8 = np.where(nu_rk8_diff < -180)[0]
    
    print(f"Periapsis passages (analytical): {len(wraparound_analytical)}")
    print(f"Periapsis passages (numerical): {len(wraparound_rk8)}")
    
    if len(wraparound_analytical) > 0:
        print("Analytical wraparound times:")
        for idx in wraparound_analytical:
            print(f"  At t={time_orbits[idx]:.3f} orbits: {nu_analytical_deg[idx]:.1f}° -> {nu_analytical_deg[idx+1]:.1f}°")
    
    # Create comprehensive plots
    fig = plt.figure(figsize=(16, 12))
    
    # Plot 1: Raw true anomaly vs time
    ax1 = plt.subplot(3, 2, 1)
    ax1.plot(time_orbits, nu_analytical_deg, 'r-', label='Analytical', linewidth=2)
    ax1.plot(time_orbits, nu_rk8_deg, 'b--', label='RK8 Numerical', alpha=0.7, linewidth=1)
    ax1.set_xlabel('Time (orbits)')
    ax1.set_ylabel('True Anomaly (degrees)')
    ax1.set_title('True Anomaly vs Time (Raw - with wraparound)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(-10, 370)
    
    # Plot 2: Unwrapped true anomaly vs time
    ax2 = plt.subplot(3, 2, 2)
    nu_analytical_unwrapped = np.unwrap(nu_analytical)
    nu_rk8_unwrapped = np.unwrap(nu_rk8)
    ax2.plot(time_orbits, np.degrees(nu_analytical_unwrapped), 'r-', label='Analytical', linewidth=2)
    ax2.plot(time_orbits, np.degrees(nu_rk8_unwrapped), 'b--', label='RK8 Numerical', alpha=0.7, linewidth=1)
    ax2.set_xlabel('Time (orbits)')
    ax2.set_ylabel('True Anomaly (degrees, unwrapped)')
    ax2.set_title('True Anomaly vs Time (Unwrapped - continuous)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Rate of change of true anomaly
    ax3 = plt.subplot(3, 2, 3)
    dt = np.diff(times_rk8)
    dnu_dt_analytical = nu_analytical_diff / dt * 3600  # degrees per hour
    dnu_dt_rk8 = nu_rk8_diff / dt * 3600
    ax3.plot(time_orbits[1:], dnu_dt_analytical, 'r-', label='Analytical dν/dt', linewidth=2)
    ax3.plot(time_orbits[1:], dnu_dt_rk8, 'b--', label='RK8 dν/dt', alpha=0.7)
    ax3.set_xlabel('Time (orbits)')
    ax3.set_ylabel('dν/dt (degrees/hour)')
    ax3.set_title('Rate of True Anomaly Change')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.set_yscale('log')
    
    # Plot 4: True anomaly error
    ax4 = plt.subplot(3, 2, 4)
    nu_error = nu_rk8 - nu_analytical
    # Handle wraparound properly
    nu_error = np.where(nu_error > np.pi, nu_error - 2*np.pi, nu_error)
    nu_error = np.where(nu_error < -np.pi, nu_error + 2*np.pi, nu_error)
    ax4.plot(time_orbits, np.degrees(nu_error)*3600, 'g-', linewidth=2)  # Convert to arcseconds
    ax4.set_xlabel('Time (orbits)')
    ax4.set_ylabel('True Anomaly Error (arcseconds)')
    ax4.set_title('RK8 vs Analytical True Anomaly Error')
    ax4.grid(True, alpha=0.3)
    
    # Plot 5: Position magnitude vs time (to see orbital shape)
    ax5 = plt.subplot(3, 2, 5)
    r_analytical = np.linalg.norm(pos_analytical, axis=1)
    r_rk8 = np.linalg.norm(pos_rk8, axis=1)
    ax5.plot(time_orbits, r_analytical, 'r-', label='Analytical', linewidth=2)
    ax5.plot(time_orbits, r_rk8, 'b--', label='RK8 Numerical', alpha=0.7)
    ax5.set_xlabel('Time (orbits)')
    ax5.set_ylabel('Orbital Radius (km)')
    ax5.set_title('Orbital Radius vs Time')
    ax5.legend()
    ax5.grid(True, alpha=0.3)
    
    # Plot 6: True anomaly vs orbital radius (phase plot)
    ax6 = plt.subplot(3, 2, 6)
    ax6.plot(nu_analytical_deg, r_analytical, 'r-', label='Analytical', linewidth=2)
    ax6.plot(nu_rk8_deg, r_rk8, 'b--', label='RK8 Numerical', alpha=0.7)
    ax6.set_xlabel('True Anomaly (degrees)')
    ax6.set_ylabel('Orbital Radius (km)')
    ax6.set_title('Radius vs True Anomaly (Shape)')
    ax6.legend()
    ax6.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('full_orbit_true_anomaly_analysis.png', dpi=300, bbox_inches='tight')
    print(f"\nDetailed plots saved as 'full_orbit_true_anomaly_analysis.png'")
    
    # Statistics
    print(f"\nStatistics:")
    print(f"Max true anomaly error: {np.degrees(np.max(np.abs(nu_error)))*3600:.3f} arcseconds")
    print(f"RMS true anomaly error: {np.degrees(np.sqrt(np.mean(nu_error**2)))*3600:.3f} arcseconds")
    
    # Check if the issue is in the plotting
    print(f"\nDiagnostics:")
    print(f"Max rate of ν change (analytical): {np.max(np.abs(dnu_dt_analytical)):.1f}°/hour")
    print(f"Min rate of ν change (analytical): {np.min(np.abs(dnu_dt_analytical)):.3f}°/hour") 
    print(f"Rate ratio (max/min): {np.max(np.abs(dnu_dt_analytical))/np.min(np.abs(dnu_dt_analytical)):.1f}:1")
    
    return time_orbits, nu_analytical, nu_rk8

if __name__ == "__main__":
    test_full_orbit_true_anomaly()