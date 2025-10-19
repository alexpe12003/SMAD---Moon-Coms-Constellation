"""
Main Simulation Script for Orbit B Numerical Integration Study
=============================================================

This script runs the complete numerical integration study for Orbit B,
comparing RK4 numerical integration with analytical solutions and
analyzing convergence properties.

Author: Generated for SMAD Moon Communications Constellation Study
Date: October 2025
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import os
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import our modules
from core.Constants import (
    get_orbit_b_parameters, 
    INTEGRATION_STEPS, 
    NUM_ORBITS_STUDY,
    RESULTS_DIR,
    print_orbit_summary
)
from core.Kepler import analytical_propagation, position_velocity_from_elements, analytical_true_anomaly_propagation
from core.RK4 import propagate_orbit_rk4, calculate_orbital_energy, calculate_angular_momentum, calculate_true_anomaly_from_state

# =============================================================================
# MAIN SIMULATION CLASS
# =============================================================================

class OrbitSimulation:
    """Main class for running orbit propagation studies."""
    
    def __init__(self):
        """Initialize the simulation with Orbit B parameters."""
        self.orbit_params = get_orbit_b_parameters()
        self.results = {}
        
        # Create results directory
        Path(RESULTS_DIR).mkdir(exist_ok=True)
        
        print("Orbit B Numerical Integration Study")
        print("=" * 50)
        print_orbit_summary()
    
    def setup_initial_conditions(self):
        """Set up initial conditions for the simulation."""
        # Get orbital parameters
        a = self.orbit_params['semimajor_axis']
        e = self.orbit_params['eccentricity']
        i = self.orbit_params['inclination']
        raan = self.orbit_params['raan']
        arg_p = self.orbit_params['arg_periapsis']
        nu_0 = self.orbit_params['true_anomaly_0']
        
        # Calculate initial position and velocity
        r0, v0 = position_velocity_from_elements(a, e, i, raan, arg_p, nu_0)
        
        # Store initial state
        self.initial_state = np.concatenate([r0, v0])
        self.initial_position = r0
        self.initial_velocity = v0
        
        print(f"\nInitial Conditions:")
        print(f"Position: [{r0[0]:>10.3f}, {r0[1]:>10.3f}, {r0[2]:>10.3f}] km")
        print(f"Velocity: [{v0[0]:>10.6f}, {v0[1]:>10.6f}, {v0[2]:>10.6f}] km/s")
        print(f"Initial orbital radius: {np.linalg.norm(r0):>10.3f} km")
        print(f"Initial speed: {np.linalg.norm(v0):>10.6f} km/s")
        
        return self.initial_state
    
    def run_convergence_study(self):
        """Run convergence study with different step sizes."""
        print(f"\nRunning convergence study...")
        print(f"Integration steps to test: {INTEGRATION_STEPS}")
        print(f"Number of orbits: {NUM_ORBITS_STUDY}")
        
        # Time span for integration
        period = self.orbit_params['orbital_period']
        t_end = NUM_ORBITS_STUDY * period
        
        # Get analytical solution as reference with 1s time step
        analytical_step_size = 1.0  # 1 second time step
        time_ref = np.arange(0, t_end + analytical_step_size, analytical_step_size)
        pos_analytical, vel_analytical = analytical_propagation(
            self.orbit_params, time_ref
        )
        
        # Get analytical true anomaly as reference
        true_anomaly_analytical = analytical_true_anomaly_propagation(
            self.orbit_params, time_ref
        )
        
        convergence_results = {}
        
        for step_size in INTEGRATION_STEPS:
            print(f"\n  Testing step size: {step_size} s")
            
            start_time = time.time()
            
            # Run RK4 integration
            times_rk4, pos_rk4, vel_rk4 = propagate_orbit_rk4(
                self.initial_state,
                (0, t_end),
                step_size
            )
            
            computation_time = time.time() - start_time
            
            # Interpolate analytical solution to RK4 time points
            pos_analytical_interp = np.array([
                np.interp(times_rk4, time_ref, pos_analytical[:, i])
                for i in range(3)
            ]).T
            
            vel_analytical_interp = np.array([
                np.interp(times_rk4, time_ref, vel_analytical[:, i])
                for i in range(3)
            ]).T
            
            # Calculate errors
            position_errors = np.linalg.norm(pos_rk4 - pos_analytical_interp, axis=1)
            velocity_errors = np.linalg.norm(vel_rk4 - vel_analytical_interp, axis=1)
            
            # Calculate energy conservation
            energy_rk4 = calculate_orbital_energy(pos_rk4, vel_rk4)
            energy_error = np.abs(energy_rk4 - energy_rk4[0])
            
            # Calculate angular momentum conservation
            h_rk4 = calculate_angular_momentum(pos_rk4, vel_rk4)
            h_error = np.abs(h_rk4 - h_rk4[0])
            
            # Calculate true anomaly from RK4 results
            true_anomaly_rk4 = calculate_true_anomaly_from_state(
                pos_rk4, vel_rk4, self.orbit_params
            )
            
            # Interpolate analytical true anomaly to RK4 time points
            true_anomaly_analytical_interp = np.interp(times_rk4, time_ref, true_anomaly_analytical)
            
            # Calculate true anomaly errors (handle 2π wraparound)
            true_anomaly_diff = true_anomaly_rk4 - true_anomaly_analytical_interp
            # Correct for wraparound: if difference > π, subtract 2π; if < -π, add 2π
            true_anomaly_diff = np.where(true_anomaly_diff > np.pi, 
                                       true_anomaly_diff - 2*np.pi, true_anomaly_diff)
            true_anomaly_diff = np.where(true_anomaly_diff < -np.pi, 
                                       true_anomaly_diff + 2*np.pi, true_anomaly_diff)
            true_anomaly_errors = np.abs(true_anomaly_diff)
            
            # Store results
            convergence_results[step_size] = {
                'times': times_rk4,
                'positions': pos_rk4,
                'velocities': vel_rk4,
                'position_errors': position_errors,
                'velocity_errors': velocity_errors,
                'energy_errors': energy_error,
                'angular_momentum_errors': h_error,
                'true_anomaly_rk4': true_anomaly_rk4,
                'true_anomaly_errors': true_anomaly_errors,
                'max_position_error': np.max(position_errors),
                'max_velocity_error': np.max(velocity_errors),
                'max_energy_error': np.max(energy_error),
                'max_h_error': np.max(h_error),
                'max_true_anomaly_error': np.max(true_anomaly_errors),
                'final_position_error': position_errors[-1],
                'final_velocity_error': velocity_errors[-1],
                'final_true_anomaly_error': true_anomaly_errors[-1],
                'computation_time': computation_time,
                'num_steps': len(times_rk4)
            }
            
            print(f"    Steps: {len(times_rk4):>6d}, Time: {computation_time:>6.3f} s")
            print(f"    Max position error:     {np.max(position_errors):>10.3e} km")
            print(f"    Max velocity error:     {np.max(velocity_errors):>10.3e} km/s")
            print(f"    Max energy error:       {np.max(energy_error):>10.3e} km²/s²")
            print(f"    Max true anomaly error: {np.max(true_anomaly_errors):>10.3e} rad ({np.degrees(np.max(true_anomaly_errors)):>8.3e}°)")
        
        self.results['convergence'] = convergence_results
        self.results['analytical_reference'] = {
            'times': time_ref,
            'positions': pos_analytical,
            'velocities': vel_analytical,
            'true_anomaly': true_anomaly_analytical
        }
        
        return convergence_results
    
    def analyze_results(self):
        """Analyze and print convergence results."""
        print(f"\n{'='*60}")
        print("CONVERGENCE ANALYSIS RESULTS")
        print(f"{'='*60}")
        
        convergence_data = self.results['convergence']
        
        # Print summary table
        print(f"{'Step Size':>10} {'Steps':>8} {'Time':>8} {'Max Pos Err':>12} {'Max Vel Err':>12} {'Energy Err':>12} {'True Anom Err':>14}")
        print(f"{'(s)':>10} {'':>8} {'(s)':>8} {'(km)':>12} {'(km/s)':>12} {'(km²/s²)':>12} {'(rad)':>14}")
        print("-" * 94)
        
        for step_size in INTEGRATION_STEPS:
            result = convergence_data[step_size]
            print(f"{step_size:>10.1f} {result['num_steps']:>8d} "
                  f"{result['computation_time']:>8.3f} "
                  f"{result['max_position_error']:>12.3e} "
                  f"{result['max_velocity_error']:>12.3e} "
                  f"{result['max_energy_error']:>12.3e} "
                  f"{result['max_true_anomaly_error']:>14.3e}")
        
        # Analyze convergence order
        step_sizes = np.array(INTEGRATION_STEPS)
        max_pos_errors = np.array([convergence_data[h]['max_position_error'] 
                                  for h in INTEGRATION_STEPS])
        
        # Fit convergence order (log-log slope)
        valid_indices = (max_pos_errors > 0) & (step_sizes > 0)
        if np.sum(valid_indices) >= 2:
            log_h = np.log(step_sizes[valid_indices])
            log_err = np.log(max_pos_errors[valid_indices])
            convergence_order = np.polyfit(log_h, log_err, 1)[0]
            
            print(f"\nEstimated convergence order: {convergence_order:.2f}")
            print(f"Expected for RK4: ~4.0")
            
    def save_results(self):
        """Save results to files."""
        print(f"\nSaving results to {RESULTS_DIR}/")
        
        convergence_data = self.results['convergence']
        
        # Save convergence summary
        with open(f"{RESULTS_DIR}/convergence_summary.txt", 'w') as f:
            f.write("Orbit B Convergence Study Results\n")
            f.write("=" * 40 + "\n\n")
            f.write("Orbital Parameters:\n")
            for key, value in self.orbit_params.items():
                f.write(f"  {key}: {value}\n")
            f.write(f"\nIntegration Steps: {INTEGRATION_STEPS}\n")
            f.write(f"Number of Orbits: {NUM_ORBITS_STUDY}\n\n")
            
            f.write(f"{'Step Size':>10} {'Steps':>8} {'Time':>8} {'Max Pos Err':>12} {'Max Vel Err':>12} {'True Anom Err':>14}\n")
            f.write(f"{'(s)':>10} {'':>8} {'(s)':>8} {'(km)':>12} {'(km/s)':>12} {'(rad)':>14}\n")
            f.write("-" * 74 + "\n")
            
            for step_size in INTEGRATION_STEPS:
                result = convergence_data[step_size]
                f.write(f"{step_size:>10.1f} {result['num_steps']:>8d} "
                       f"{result['computation_time']:>8.3f} "
                       f"{result['max_position_error']:>12.3e} "
                       f"{result['max_velocity_error']:>12.3e} "
                       f"{result['max_true_anomaly_error']:>14.3e}\n")
        
        # Save detailed results for each step size
        for step_size in INTEGRATION_STEPS:
            result = convergence_data[step_size]
            filename = f"{RESULTS_DIR}/rk4_results_step_{step_size:g}s.npz"
            
            np.savez(filename,
                     times=result['times'],
                     positions=result['positions'],
                     velocities=result['velocities'],
                     position_errors=result['position_errors'],
                     velocity_errors=result['velocity_errors'],
                     energy_errors=result['energy_errors'],
                     true_anomaly_rk4=result['true_anomaly_rk4'],
                     true_anomaly_errors=result['true_anomaly_errors'],
                     step_size=step_size)
        
        print(f"Results saved successfully!")

# =============================================================================
# PLOTTING UTILITIES
# =============================================================================

def create_plots(simulation):
    """Create comprehensive plots of the simulation results."""
    print(f"\nCreating plots...")
    
    convergence_data = simulation.results['convergence']
    analytical_data = simulation.results['analytical_reference']
    
    # Set up plotting style
    plt.style.use('default')
    
    # 1. Create comprehensive figure with 2x3 subplot layout
    fig = plt.figure(figsize=(20, 12))
    
    # 3D trajectory plot
    ax = fig.add_subplot(2, 3, 1, projection='3d')
    
    # Plot analytical solution
    pos_analytical = analytical_data['positions']
    ax.plot(pos_analytical[:, 0], pos_analytical[:, 1], pos_analytical[:, 2],
            'k-', label='Analytical', linewidth=2)
    
    # Plot RK4 solutions for different step sizes
    colors = plt.cm.viridis(np.linspace(0, 1, len(INTEGRATION_STEPS)))
    for i, step_size in enumerate(INTEGRATION_STEPS[::2]):  # Every other step for clarity
        result = convergence_data[step_size]
        pos_rk4 = result['positions']
        ax.plot(pos_rk4[:, 0], pos_rk4[:, 1], pos_rk4[:, 2],
                '--', color=colors[i*2], label=f'RK4 (Δt={step_size}s)', alpha=0.7)
    
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.set_title('Orbital Trajectories')
    ax.legend()
    
    # 2. Convergence plot
    ax = fig.add_subplot(2, 3, 2)
    step_sizes = np.array(INTEGRATION_STEPS)
    max_pos_errors = np.array([convergence_data[h]['max_position_error'] 
                              for h in INTEGRATION_STEPS])
    
    ax.loglog(step_sizes, max_pos_errors, 'bo-', label='Position Error')
    ax.loglog(step_sizes, step_sizes**4 * max_pos_errors[0] / step_sizes[0]**4, 
              'r--', label='4th Order Reference')
    ax.set_xlabel('Step Size (s)')
    ax.set_ylabel('Maximum Position Error (km)')
    ax.set_title('Convergence Study')
    ax.legend()
    ax.grid(True)
    
    # 3. True Anomaly Error over Time
    ax = fig.add_subplot(2, 3, 3)
    for i, step_size in enumerate(INTEGRATION_STEPS[::2]):
        result = convergence_data[step_size]
        times = result['times'] / simulation.orbit_params['orbital_period']
        true_anomaly_errors = result['true_anomaly_errors']
        ax.semilogy(times, true_anomaly_errors, 
                   color=colors[i*2], label=f'Δt={step_size}s')
    
    ax.set_xlabel('Time (orbital periods)')
    ax.set_ylabel('True Anomaly Error (rad)')
    ax.set_title('True Anomaly Error Evolution')
    ax.legend()
    ax.grid(True)
    
    # 4. True Anomaly Convergence
    ax = fig.add_subplot(2, 3, 4)
    max_true_anomaly_errors = np.array([convergence_data[h]['max_true_anomaly_error'] 
                                      for h in INTEGRATION_STEPS])
    
    ax.loglog(step_sizes, max_true_anomaly_errors, 'mo-', label='True Anomaly Error')
    ax.loglog(step_sizes, step_sizes**4 * max_true_anomaly_errors[0] / step_sizes[0]**4, 
              'r--', label='4th Order Reference')
    ax.set_xlabel('Step Size (s)')
    ax.set_ylabel('Maximum True Anomaly Error (rad)')
    ax.set_title('True Anomaly Convergence')
    ax.legend()
    ax.grid(True)
    
    # 5. Energy conservation
    ax = fig.add_subplot(2, 3, 5)
    for i, step_size in enumerate(INTEGRATION_STEPS[::2]):
        result = convergence_data[step_size]
        times = result['times'] / simulation.orbit_params['orbital_period']
        energy_errors = result['energy_errors']
        ax.semilogy(times, np.abs(energy_errors), 
                   color=colors[i*2], label=f'Δt={step_size}s')
    
    ax.set_xlabel('Time (orbital periods)')
    ax.set_ylabel('|Energy Error| (km²/s²)')
    ax.set_title('Energy Conservation')
    ax.legend()
    ax.grid(True)
    
    # 6. Computation efficiency
    ax = fig.add_subplot(2, 3, 6)
    comp_times = np.array([convergence_data[h]['computation_time'] 
                          for h in INTEGRATION_STEPS])
    
    ax.loglog(step_sizes, comp_times, 'go-', label='Computation Time')
    ax.loglog(step_sizes, max_pos_errors, 'bo-', label='Position Error')
    ax.set_xlabel('Step Size (s)')
    ax.set_ylabel('Time (s) / Error (km)')
    ax.set_title('Efficiency Analysis')
    ax.legend()
    ax.grid(True)
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/orbit_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Create dedicated true anomaly analysis plot
    create_true_anomaly_plot(simulation)
    
    print(f"Plots saved to {RESULTS_DIR}/orbit_analysis.png and true_anomaly_analysis.png")

def create_true_anomaly_plot(simulation):
    """Create dedicated plot for true anomaly analysis."""
    print(f"Creating detailed true anomaly analysis plot...")
    
    convergence_data = simulation.results['convergence']
    analytical_data = simulation.results['analytical_reference']
    
    # Create figure for true anomaly analysis
    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(INTEGRATION_STEPS)))
    
    # 1. True anomaly comparison (analytical vs numerical)
    ax = axes[0, 0]
    
    # Plot analytical true anomaly
    times_analytical = analytical_data['times'] / simulation.orbit_params['orbital_period']
    true_anomaly_analytical = analytical_data['true_anomaly']
    ax.plot(times_analytical, np.degrees(true_anomaly_analytical), 
            'k-', label='Analytical', linewidth=2)
    
    # Plot RK4 results for selected step sizes
    for i, step_size in enumerate(INTEGRATION_STEPS[::3]):  # Every third for clarity
        result = convergence_data[step_size]
        times = result['times'] / simulation.orbit_params['orbital_period']
        true_anomaly_rk4 = result['true_anomaly_rk4']
        ax.plot(times, np.degrees(true_anomaly_rk4), 
                '--', color=colors[i*3], label=f'RK4 (Δt={step_size}s)', alpha=0.8)
    
    ax.set_xlabel('Time (orbital periods)')
    ax.set_ylabel('True Anomaly (degrees)')
    ax.set_title('True Anomaly: Analytical vs RK4')
    ax.legend()
    ax.grid(True)
    
    # 2. True anomaly error vs time for all step sizes
    ax = axes[0, 1]
    for i, step_size in enumerate(INTEGRATION_STEPS):
        result = convergence_data[step_size]
        times = result['times'] / simulation.orbit_params['orbital_period']
        true_anomaly_errors = result['true_anomaly_errors']
        ax.semilogy(times, np.degrees(true_anomaly_errors), 
                   color=colors[i], label=f'Δt={step_size}s', alpha=0.7)
    
    ax.set_xlabel('Time (orbital periods)')
    ax.set_ylabel('True Anomaly Error (degrees)')
    ax.set_title('True Anomaly Error Evolution')
    ax.legend()
    ax.grid(True)
    
    # 3. Convergence analysis for true anomaly
    ax = axes[1, 0]
    step_sizes = np.array(INTEGRATION_STEPS)
    max_true_anomaly_errors = np.array([convergence_data[h]['max_true_anomaly_error'] 
                                      for h in INTEGRATION_STEPS])
    final_true_anomaly_errors = np.array([convergence_data[h]['final_true_anomaly_error'] 
                                        for h in INTEGRATION_STEPS])
    
    ax.loglog(step_sizes, max_true_anomaly_errors, 'mo-', label='Maximum Error')
    ax.loglog(step_sizes, final_true_anomaly_errors, 'co-', label='Final Error')
    ax.loglog(step_sizes, step_sizes**4 * max_true_anomaly_errors[0] / step_sizes[0]**4, 
              'r--', label='4th Order Reference')
    ax.set_xlabel('Step Size (s)')
    ax.set_ylabel('True Anomaly Error (rad)')
    ax.set_title('True Anomaly Convergence Study')
    ax.legend()
    ax.grid(True)
    
    # 4. Error statistics comparison
    ax = axes[1, 1]
    x_pos = np.arange(len(INTEGRATION_STEPS))
    width = 0.35
    
    max_pos_errors = np.array([convergence_data[h]['max_position_error'] 
                              for h in INTEGRATION_STEPS])
    
    # Normalize errors for comparison (to show relative behavior)
    pos_err_norm = max_pos_errors / np.max(max_pos_errors)
    anom_err_norm = max_true_anomaly_errors / np.max(max_true_anomaly_errors)
    
    ax.bar(x_pos - width/2, pos_err_norm, width, label='Position Error (norm.)', alpha=0.8)
    ax.bar(x_pos + width/2, anom_err_norm, width, label='True Anomaly Error (norm.)', alpha=0.8)
    
    ax.set_xlabel('Integration Step Size')
    ax.set_ylabel('Normalized Error')
    ax.set_title('Position vs True Anomaly Error Comparison')
    ax.set_xticks(x_pos)
    ax.set_xticklabels([f'{h}s' for h in INTEGRATION_STEPS])
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/true_anomaly_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()

# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main():
    """Run the complete orbit simulation study."""
    print("Starting Orbit B Numerical Integration Study...")
    print(f"Current working directory: {os.getcwd()}")
    
    # Create simulation object
    sim = OrbitSimulation()
    
    # Set up initial conditions
    sim.setup_initial_conditions()
    
    # Run convergence study
    sim.run_convergence_study()
    
    # Analyze results
    sim.analyze_results()
    
    # Save results
    sim.save_results()
    
    # Create plots
    create_plots(sim)
    
    print("\nSimulation completed successfully!")
    print(f"Results saved in '{RESULTS_DIR}/' directory")

if __name__ == "__main__":
    main()