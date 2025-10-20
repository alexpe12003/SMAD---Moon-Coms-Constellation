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
from core.Kepler import (
    analytical_propagation, 
    position_velocity_from_elements, 
    analytical_true_anomaly_propagation
)
from core.RK4 import propagate_orbit_rk4, calculate_orbital_energy, calculate_angular_momentum, calculate_true_anomaly_from_state
from core.RK8 import propagate_orbit_rk8
from plotting_utils import create_all_plots

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
        print(f"Running convergence study...")
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
        
        # Test both RK4 and RK8 methods
        integration_methods = {
            'RK4': propagate_orbit_rk4,
            'RK8': propagate_orbit_rk8
        }
        
        for method_name, integration_func in integration_methods.items():
            print(f"\n--- Testing {method_name} Method ---")
            convergence_results[method_name] = {}
            
            for step_size in INTEGRATION_STEPS:
                print(f"\n  {method_name} - Testing step size: {step_size} s")
                
                start_time = time.time()
                
                # Run integration
                times_method, pos_method, vel_method = integration_func(
                    self.initial_state,
                    (0, t_end),
                    step_size
                )
                
                computation_time = time.time() - start_time
                
                # Interpolate analytical solution to method time points
                pos_analytical_interp = np.array([
                    np.interp(times_method, time_ref, pos_analytical[:, i])
                    for i in range(3)
                ]).T
                
                vel_analytical_interp = np.array([
                    np.interp(times_method, time_ref, vel_analytical[:, i])
                    for i in range(3)
                ]).T
                # Calculate errors
                position_errors = np.linalg.norm(pos_method - pos_analytical_interp, axis=1)
                velocity_errors = np.linalg.norm(vel_method - vel_analytical_interp, axis=1)
                
                # Calculate energy conservation
                energy_method = calculate_orbital_energy(pos_method, vel_method)
                energy_error = np.abs(energy_method - energy_method[0])
                
                # Calculate angular momentum conservation
                h_method = calculate_angular_momentum(pos_method, vel_method)
                h_error = np.abs(h_method - h_method[0])
                
                # Calculate true anomaly from method results
                true_anomaly_method = calculate_true_anomaly_from_state(
                    pos_method, vel_method, self.orbit_params
                )
                
                # Interpolate analytical true anomaly to method time points
                true_anomaly_analytical_interp = np.interp(times_method, time_ref, true_anomaly_analytical)
                
                # Calculate true anomaly errors (handle 2π wraparound)
                true_anomaly_diff = true_anomaly_method - true_anomaly_analytical_interp
                # Correct for wraparound: if difference > π, subtract 2π; if < -π, add 2π
                true_anomaly_diff = np.where(true_anomaly_diff > np.pi, 
                                           true_anomaly_diff - 2*np.pi, true_anomaly_diff)
                true_anomaly_diff = np.where(true_anomaly_diff < -np.pi, 
                                           true_anomaly_diff + 2*np.pi, true_anomaly_diff)
                true_anomaly_errors = np.abs(true_anomaly_diff)
                
                # Store results
                convergence_results[method_name][step_size] = {
                    'times': times_method,
                    'positions': pos_method,
                    'velocities': vel_method,
                    'position_errors': position_errors,
                    'velocity_errors': velocity_errors,
                    'energy_errors': energy_error,
                    'angular_momentum_errors': h_error,
                    'true_anomaly_method': true_anomaly_method,
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
                    'num_steps': len(times_method)
                }
                
                print(f"    Steps: {len(times_method):>6d}, Time: {computation_time:>6.3f} s")
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
        print(f"\n{'='*80}")
        print("CONVERGENCE ANALYSIS RESULTS")
        print(f"{'='*80}")
        
        convergence_data = self.results['convergence']
        analytical_data = self.results['analytical_reference']
        
        print(f"\nAnalytical Reference Solution:")
        print(f"  Time points: {len(analytical_data['times']):,}")
        
        for method_name in ['RK4', 'RK8']:
            print(f"\n{method_name} Method Results:")
            print(f"{'='*40}")
            
            # Print summary table for this method
            print(f"{'Step Size':>10} {'Steps':>8} {'Time':>8} {'Max Pos Err':>12} {'Max Vel Err':>12} {'Energy Err':>12} {'True Anom Err':>14}")
            print(f"{'(s)':>10} {'':>8} {'(s)':>8} {'(km)':>12} {'(km/s)':>12} {'(km²/s²)':>12} {'(rad)':>14}")
            print("-" * 94)
            
            for step_size in INTEGRATION_STEPS:
                result = convergence_data[method_name][step_size]
                print(f"{step_size:>10.1f} {result['num_steps']:>8d} "
                      f"{result['computation_time']:>8.3f} "
                      f"{result['max_position_error']:>12.3e} "
                      f"{result['max_velocity_error']:>12.3e} "
                      f"{result['max_energy_error']:>12.3e} "
                      f"{result['max_true_anomaly_error']:>14.3e}")
            
            # Analyze convergence order for this method
            step_sizes = np.array(INTEGRATION_STEPS)
            max_pos_errors = np.array([convergence_data[method_name][h]['max_position_error'] 
                                      for h in INTEGRATION_STEPS])
            
            # Fit convergence order (log-log slope)
            valid_indices = (max_pos_errors > 0) & (step_sizes > 0)
            if np.sum(valid_indices) >= 2:
                log_h = np.log(step_sizes[valid_indices])
                log_err = np.log(max_pos_errors[valid_indices])
                convergence_order = np.polyfit(log_h, log_err, 1)[0]
                
                print(f"\nEstimated convergence order: {convergence_order:.2f}")
                if method_name == 'RK4':
                    print(f"Expected for RK4: ~4.0")
                elif method_name == 'RK8':
                    print(f"Expected for RK8: ~8.0")
            
    def save_results(self):
        """Save results to files."""
        print(f"\nSaving results to {RESULTS_DIR}/")
        
        convergence_data = self.results['convergence']
        
        # Save convergence summary for both methods
        with open(f"{RESULTS_DIR}/convergence_summary.txt", 'w') as f:
            f.write("Orbit B Convergence Study Results\n")
            f.write("=" * 40 + "\n\n")
            f.write("Orbital Parameters:\n")
            for key, value in self.orbit_params.items():
                f.write(f"  {key}: {value}\n")
            f.write(f"\nIntegration Steps: {INTEGRATION_STEPS}\n")
            f.write(f"Number of Orbits: {NUM_ORBITS_STUDY}\n\n")
            
            for method_name in ['RK4', 'RK8']:
                f.write(f"\n{method_name} Method Results:\n")
                f.write("=" * 30 + "\n")
                f.write(f"{'Step Size':>10} {'Steps':>8} {'Time':>8} {'Max Pos Err':>12} {'Max Vel Err':>12} {'True Anom Err':>14}\n")
                f.write(f"{'(s)':>10} {'':>8} {'(s)':>8} {'(km)':>12} {'(km/s)':>12} {'(rad)':>14}\n")
                f.write("-" * 74 + "\n")
                
                for step_size in INTEGRATION_STEPS:
                    result = convergence_data[method_name][step_size]
                    f.write(f"{step_size:>10.1f} {result['num_steps']:>8d} "
                           f"{result['computation_time']:>8.3f} "
                           f"{result['max_position_error']:>12.3e} "
                           f"{result['max_velocity_error']:>12.3e} "
                           f"{result['max_true_anomaly_error']:>14.3e}\n")
        
        # Save detailed results for each method and step size
        for method_name in ['RK4', 'RK8']:
            for step_size in INTEGRATION_STEPS:
                result = convergence_data[method_name][step_size]
                filename = f"{RESULTS_DIR}/{method_name.lower()}_results_step_{step_size:g}s.npz"
                
                np.savez(filename,
                         times=result['times'],
                         positions=result['positions'],
                         velocities=result['velocities'],
                         position_errors=result['position_errors'],
                         velocity_errors=result['velocity_errors'],
                         energy_errors=result['energy_errors'],
                         true_anomaly_method=result['true_anomaly_method'],
                         true_anomaly_errors=result['true_anomaly_errors'],
                         step_size=step_size,
                         method=method_name)
        
        print(f"Results saved successfully!")

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
    create_all_plots(sim)
    
    print("\nSimulation completed successfully!")
    print(f"Results saved in '{RESULTS_DIR}/' directory")

if __name__ == "__main__":
    main()