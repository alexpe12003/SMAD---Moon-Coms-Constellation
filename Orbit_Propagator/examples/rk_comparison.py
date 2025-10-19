"""
RK4 vs RK8 Comparison Script
============================

This script demonstrates the difference between RK4 and RK8 integration
methods for the same orbital mechanics problem.

Author: Generated for SMAD Moon Communications Constellation Study
Date: October 2025
"""

import numpy as np
import matplotlib.pyplot as plt
import time
from Constants import get_orbit_b_parameters
from Kepler import position_velocity_from_elements, analytical_propagation
from RK4 import propagate_orbit_rk4
from RK8 import propagate_orbit_rk8, compare_rk4_rk8

def run_rk_comparison():
    """Run a comparison between RK4 and RK8 methods."""
    
    print("RK4 vs RK8 Comparison Study")
    print("=" * 40)
    
    # Get orbital parameters
    orbit_params = get_orbit_b_parameters()
    
    # Set up initial conditions at periapsis
    a = orbit_params['semimajor_axis']
    e = orbit_params['eccentricity']
    i = orbit_params['inclination']
    raan = orbit_params['raan']
    arg_p = orbit_params['arg_periapsis']
    nu_0 = 0.0  # Start at periapsis
    
    # Calculate initial state
    r0, v0 = position_velocity_from_elements(a, e, i, raan, arg_p, nu_0)
    initial_state = np.concatenate([r0, v0])
    
    print(f"Initial position: [{r0[0]:>10.3f}, {r0[1]:>10.3f}, {r0[2]:>10.3f}] km")
    print(f"Initial velocity: [{v0[0]:>10.6f}, {v0[1]:>10.6f}, {v0[2]:>10.6f}] km/s")
    
    # Integration parameters
    period = orbit_params['orbital_period']
    t_end = period  # One orbit
    step_size = 60.0  # 1 minute steps
    
    print(f"\nIntegration setup:")
    print(f"Time span: 0 to {t_end/3600:.2f} hours (1 orbit)")
    print(f"Step size: {step_size} s")
    
    # Run comparison
    print(f"\nRunning integration...")
    
    # RK4 integration
    start_time = time.time()
    times_rk4, pos_rk4, vel_rk4 = propagate_orbit_rk4(
        initial_state, (0, t_end), step_size
    )
    rk4_time = time.time() - start_time
    
    # RK8 integration
    start_time = time.time()
    times_rk8, pos_rk8, vel_rk8 = propagate_orbit_rk8(
        initial_state, (0, t_end), step_size
    )
    rk8_time = time.time() - start_time
    
    # Analytical solution
    time_analytical = np.linspace(0, t_end, 1000)
    pos_analytical, vel_analytical = analytical_propagation(
        {**orbit_params, 'true_anomaly_0': nu_0}, time_analytical
    )
    
    # Interpolate analytical solution to integration time points
    pos_analytical_rk4 = np.array([
        np.interp(times_rk4, time_analytical, pos_analytical[:, i])
        for i in range(3)
    ]).T
    
    pos_analytical_rk8 = np.array([
        np.interp(times_rk8, time_analytical, pos_analytical[:, i])
        for i in range(3)
    ]).T
    
    # Calculate errors
    pos_error_rk4 = np.linalg.norm(pos_rk4 - pos_analytical_rk4, axis=1)
    pos_error_rk8 = np.linalg.norm(pos_rk8 - pos_analytical_rk8, axis=1)
    
    # Print results
    print(f"\nResults:")
    print(f"{'Method':<8} {'Steps':<8} {'Time (s)':<10} {'Max Error (km)':<15} {'Final Error (km)':<18}")
    print("-" * 65)
    print(f"{'RK4':<8} {len(times_rk4):<8} {rk4_time:<10.3f} {np.max(pos_error_rk4):<15.3e} {pos_error_rk4[-1]:<18.3e}")
    print(f"{'RK8':<8} {len(times_rk8):<8} {rk8_time:<10.3f} {np.max(pos_error_rk8):<15.3e} {pos_error_rk8[-1]:<18.3e}")
    
    # Error improvement
    error_improvement = np.max(pos_error_rk4) / np.max(pos_error_rk8)
    time_ratio = rk8_time / rk4_time
    
    print(f"\nError improvement (RK8 vs RK4): {error_improvement:.1f}x better")
    print(f"Computation time ratio (RK8/RK4): {time_ratio:.1f}x slower")
    print(f"Efficiency (error improvement per time cost): {error_improvement/time_ratio:.1f}")
    
    # Create comparison plot
    create_comparison_plot(times_rk4, pos_error_rk4, times_rk8, pos_error_rk8, 
                          pos_rk4, pos_rk8, pos_analytical_rk4, period)
    
    return {
        'rk4_error': np.max(pos_error_rk4),
        'rk8_error': np.max(pos_error_rk8),
        'rk4_time': rk4_time,
        'rk8_time': rk8_time,
        'error_improvement': error_improvement,
        'time_ratio': time_ratio
    }

def create_comparison_plot(times_rk4, error_rk4, times_rk8, error_rk8, 
                          pos_rk4, pos_rk8, pos_analytical, period):
    """Create comparison plots."""
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Convert times to orbital periods for better readability
    times_rk4_norm = times_rk4 / period
    times_rk8_norm = times_rk8 / period
    
    # 1. Error comparison over time
    axes[0, 0].semilogy(times_rk4_norm, error_rk4, 'b-', linewidth=2, label='RK4', alpha=0.8)
    axes[0, 0].semilogy(times_rk8_norm, error_rk8, 'r-', linewidth=2, label='RK8', alpha=0.8)
    axes[0, 0].set_xlabel('Time (orbital periods)')
    axes[0, 0].set_ylabel('Position Error (km)')
    axes[0, 0].set_title('Position Error Evolution')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # 2. 3D trajectory comparison
    ax = fig.add_subplot(2, 2, 2, projection='3d')
    ax.plot(pos_analytical[:, 0], pos_analytical[:, 1], pos_analytical[:, 2], 
            'k-', linewidth=3, label='Analytical', alpha=0.7)
    ax.plot(pos_rk4[:, 0], pos_rk4[:, 1], pos_rk4[:, 2], 
            'b--', linewidth=2, label='RK4', alpha=0.8)
    ax.plot(pos_rk8[:, 0], pos_rk8[:, 1], pos_rk8[:, 2], 
            'r:', linewidth=2, label='RK8', alpha=0.8)
    
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.set_title('Orbital Trajectories')
    ax.legend()
    
    # 3. Error ratio over time
    # Interpolate RK8 error to RK4 time points for comparison
    error_rk8_interp = np.interp(times_rk4, times_rk8, error_rk8)
    error_ratio = error_rk4 / (error_rk8_interp + 1e-16)  # Add small value to avoid division by zero
    
    axes[1, 0].plot(times_rk4_norm, error_ratio, 'g-', linewidth=2)
    axes[1, 0].set_xlabel('Time (orbital periods)')
    axes[1, 0].set_ylabel('Error Ratio (RK4/RK8)')
    axes[1, 0].set_title('RK4/RK8 Error Ratio')
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].set_yscale('log')
    
    # 4. Final error comparison bar chart
    methods = ['RK4', 'RK8']
    final_errors = [error_rk4[-1], error_rk8[-1]]
    
    bars = axes[1, 1].bar(methods, final_errors, color=['blue', 'red'], alpha=0.7)
    axes[1, 1].set_ylabel('Final Position Error (km)')
    axes[1, 1].set_title('Final Error Comparison')
    axes[1, 1].set_yscale('log')
    
    # Add value labels on bars
    for bar, error in zip(bars, final_errors):
        height = bar.get_height()
        axes[1, 1].text(bar.get_x() + bar.get_width()/2., height,
                        f'{error:.2e}', ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig('results/rk4_vs_rk8_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"\nComparison plot saved to results/rk4_vs_rk8_comparison.png")

def step_size_sensitivity_study():
    """Study sensitivity to step size for both methods."""
    
    print(f"\n" + "=" * 50)
    print("STEP SIZE SENSITIVITY STUDY")
    print("=" * 50)
    
    # Get orbital parameters
    orbit_params = get_orbit_b_parameters()
    
    # Initial conditions at periapsis
    r0, v0 = position_velocity_from_elements(
        orbit_params['semimajor_axis'],
        orbit_params['eccentricity'],
        orbit_params['inclination'],
        orbit_params['raan'],
        orbit_params['arg_periapsis'],
        0.0
    )
    initial_state = np.concatenate([r0, v0])
    
    # Test different step sizes
    step_sizes = [10.0, 30.0, 60.0, 120.0, 300.0]  # seconds
    period = orbit_params['orbital_period']
    t_end = period  # One orbit
    
    print(f"Testing step sizes: {step_sizes} seconds")
    print(f"Integration time: {t_end/3600:.2f} hours (1 orbit)")
    
    # Analytical reference
    time_ref = np.linspace(0, t_end, 2000)
    pos_ref, _ = analytical_propagation(
        {**orbit_params, 'true_anomaly_0': 0.0}, time_ref
    )
    
    results = {'step_sizes': step_sizes, 'rk4_errors': [], 'rk8_errors': [], 
               'rk4_times': [], 'rk8_times': []}
    
    for step_size in step_sizes:
        print(f"\nStep size: {step_size} s")
        
        # RK4
        start_time = time.time()
        times_rk4, pos_rk4, _ = propagate_orbit_rk4(initial_state, (0, t_end), step_size)
        rk4_time = time.time() - start_time
        
        # RK8
        start_time = time.time()
        times_rk8, pos_rk8, _ = propagate_orbit_rk8(initial_state, (0, t_end), step_size)
        rk8_time = time.time() - start_time
        
        # Interpolate reference to integration points
        pos_ref_rk4 = np.array([
            np.interp(times_rk4, time_ref, pos_ref[:, i]) for i in range(3)
        ]).T
        
        pos_ref_rk8 = np.array([
            np.interp(times_rk8, time_ref, pos_ref[:, i]) for i in range(3)
        ]).T
        
        # Calculate errors
        error_rk4 = np.max(np.linalg.norm(pos_rk4 - pos_ref_rk4, axis=1))
        error_rk8 = np.max(np.linalg.norm(pos_rk8 - pos_ref_rk8, axis=1))
        
        results['rk4_errors'].append(error_rk4)
        results['rk8_errors'].append(error_rk8)
        results['rk4_times'].append(rk4_time)
        results['rk8_times'].append(rk8_time)
        
        print(f"  RK4: Error = {error_rk4:.3e} km, Time = {rk4_time:.3f} s")
        print(f"  RK8: Error = {error_rk8:.3e} km, Time = {rk8_time:.3f} s")
        print(f"  Improvement: {error_rk4/error_rk8:.1f}x better accuracy")
    
    # Plot results
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Error vs step size
    axes[0].loglog(step_sizes, results['rk4_errors'], 'bo-', linewidth=2, markersize=8, label='RK4')
    axes[0].loglog(step_sizes, results['rk8_errors'], 'ro-', linewidth=2, markersize=8, label='RK8')
    axes[0].set_xlabel('Step Size (s)')
    axes[0].set_ylabel('Maximum Position Error (km)')
    axes[0].set_title('Convergence Study')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # Efficiency plot (error vs computation time)
    axes[1].loglog(results['rk4_times'], results['rk4_errors'], 'bo-', linewidth=2, markersize=8, label='RK4')
    axes[1].loglog(results['rk8_times'], results['rk8_errors'], 'ro-', linewidth=2, markersize=8, label='RK8')
    axes[1].set_xlabel('Computation Time (s)')
    axes[1].set_ylabel('Maximum Position Error (km)')
    axes[1].set_title('Accuracy vs Efficiency')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('results/step_size_sensitivity.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"\nStep size sensitivity plot saved to results/step_size_sensitivity.png")
    
    return results

if __name__ == "__main__":
    # Ensure results directory exists
    import os
    os.makedirs('results', exist_ok=True)
    
    # Run basic comparison
    basic_results = run_rk_comparison()
    
    # Run step size sensitivity study
    sensitivity_results = step_size_sensitivity_study()
    
    print(f"\n" + "=" * 50)
    print("SUMMARY")
    print("=" * 50)
    print(f"RK8 provides {basic_results['error_improvement']:.1f}x better accuracy")
    print(f"at {basic_results['time_ratio']:.1f}x computational cost")
    print(f"Overall efficiency gain: {basic_results['error_improvement']/basic_results['time_ratio']:.1f}x")
    print("\nFor high-precision applications, RK8 is recommended.")
    print("For real-time applications, RK4 may be sufficient.")