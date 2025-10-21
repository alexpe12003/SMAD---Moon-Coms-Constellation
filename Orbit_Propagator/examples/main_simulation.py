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
    print_orbit_summary,
    MU_EARTH
)
from core.Kepler import (
    analytical_propagation, 
    position_velocity_from_elements, 
    analytical_true_anomaly_propagation
)
from core.RK4 import propagate_orbit_rk4, calculate_orbital_energy, calculate_angular_momentum, calculate_true_anomaly_from_state
from core.RK8 import propagate_orbit_rk8
from plotting_utils import create_all_plots, create_special_case_plots, create_special_case_plots_improved, create_single_method_plots

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
                    pos_method, vel_method, MU_EARTH
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
        
        # Special case: RK8 with 100s step for 10 orbits
        print("\n" + "="*60)
        print("SPECIAL CASE: RK8 with 100s step size for 10 orbits")
        print("="*60)
        
        self.run_special_case_100_orbits()
        
        return convergence_results
    
    def run_special_case_100_orbits(self):
        """Run special case: RK8 with 100s step size for 10 orbits."""
        print(f"Running RK8 integration for 10 orbits with 100s step size...")
        
        # Parameters for special case
        num_orbits_special = 10
        step_size_special = 100.0  # seconds
        period = self.orbit_params['orbital_period']
        t_end_special = num_orbits_special * period
        
        print(f"  Total simulation time: {t_end_special:.1f} s ({t_end_special/3600:.2f} hours)")
        print(f"  Expected number of steps: ~{int(t_end_special/step_size_special):,}")
        
        # Run RK8 integration for special case
        start_time = time.time()
        times_special, pos_special, vel_special = propagate_orbit_rk8(
            self.initial_state,
            (0, t_end_special),
            step_size_special
        )
        computation_time_special = time.time() - start_time
        
        print(f"  Actual steps completed: {len(times_special):,}")
        print(f"  Computation time: {computation_time_special:.3f} s")
        
        # Get analytical solution for the same time points
        print(f"  Calculating analytical reference for comparison...")
        analytical_start = time.time()
        pos_analytical_special, vel_analytical_special = analytical_propagation(
            self.orbit_params, times_special
        )
        analytical_time_special = time.time() - analytical_start
        print(f"  Analytical computation time: {analytical_time_special:.3f} s")
        
        # Calculate errors
        position_errors_special = np.linalg.norm(pos_special - pos_analytical_special, axis=1)
        velocity_errors_special = np.linalg.norm(vel_special - vel_analytical_special, axis=1)
        
        # Calculate energy conservation
        energy_special = calculate_orbital_energy(pos_special, vel_special)
        energy_error_special = np.abs(energy_special - energy_special[0])
        
        # Calculate angular momentum conservation
        h_special = calculate_angular_momentum(pos_special, vel_special)
        h_error_special = np.abs(h_special - h_special[0])
        
        # Store special case results
        self.results['special_case_100_orbits'] = {
            'times': times_special,
            'positions': pos_special,
            'velocities': vel_special,
            'analytical_positions': pos_analytical_special,
            'analytical_velocities': vel_analytical_special,
            'position_errors': position_errors_special,
            'velocity_errors': velocity_errors_special,
            'energy_errors': energy_error_special,
            'angular_momentum_errors': h_error_special,
            'max_position_error': np.max(position_errors_special),
            'max_velocity_error': np.max(velocity_errors_special),
            'max_energy_error': np.max(energy_error_special),
            'max_h_error': np.max(h_error_special),
            'final_position_error': position_errors_special[-1],
            'final_velocity_error': velocity_errors_special[-1],
            'computation_time': computation_time_special,
            'analytical_computation_time': analytical_time_special,
            'num_steps': len(times_special),
            'step_size': step_size_special,
            'num_orbits': num_orbits_special
        }
        
        # Print summary results
        print(f"\n  Special Case Results Summary:")
        print(f"    Max position error:     {np.max(position_errors_special):>10.3e} km")
        print(f"    Final position error:   {position_errors_special[-1]:>10.3e} km")
        print(f"    Max velocity error:     {np.max(velocity_errors_special):>10.3e} km/s")
        print(f"    Max energy error:       {np.max(energy_error_special):>10.3e} km²/s²")
        print(f"    RK8 computation time:   {computation_time_special:>10.3f} s")
        print(f"    Analytical comp. time:  {analytical_time_special:>10.3f} s")
    
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
        
        # Print special case results
        if 'special_case_100_orbits' in self.results:
            special_data = self.results['special_case_100_orbits']
            print(f"\n{'='*60}")
            print("SPECIAL CASE RESULTS: RK8 - 100s step for 10 orbits")
            print(f"{'='*60}")
            print(f"Steps: {special_data['num_steps']:>8d}")
            print(f"Computation time: {special_data['computation_time']:>8.3f} s")
            print(f"Max position error: {special_data['max_position_error']:>12.3e} km")
            print(f"Final position error: {special_data['final_position_error']:>12.3e} km")
            print(f"Max velocity error: {special_data['max_velocity_error']:>12.3e} km/s")
            print(f"Max energy error: {special_data['max_energy_error']:>12.3e} km²/s²")
            
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
        
        # Save special case results
        if 'special_case_100_orbits' in self.results:
            special_data = self.results['special_case_100_orbits']
            filename = f"{RESULTS_DIR}/rk8_special_case_100orbits_100s.npz"
            
            np.savez(filename,
                     times=special_data['times'],
                     positions=special_data['positions'],
                     velocities=special_data['velocities'],
                     analytical_positions=special_data['analytical_positions'],
                     analytical_velocities=special_data['analytical_velocities'],
                     position_errors=special_data['position_errors'],
                     velocity_errors=special_data['velocity_errors'],
                     energy_errors=special_data['energy_errors'],
                     step_size=special_data['step_size'],
                     num_orbits=special_data['num_orbits'],
                     method='RK8_special_100_orbits')
        
        print(f"Results saved successfully!")

    def run_single_method_convergence_study(self, method_name):
        """Run convergence study for a single method (RK4 or RK8)."""
        print(f"Running {method_name} convergence study...")
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
        
        # Select integration method
        if method_name == 'RK4':
            integration_func = propagate_orbit_rk4
        elif method_name == 'RK8':
            integration_func = propagate_orbit_rk8
        else:
            raise ValueError(f"Unknown method: {method_name}")
        
        method_results = {}
        
        print(f"\n--- Testing {method_name} Method ---")
        
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
                pos_method, vel_method, MU_EARTH
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
            method_results[step_size] = {
                'times': times_method,
                'positions': pos_method,
                'velocities': vel_method,
                'analytical_positions': pos_analytical_interp,
                'analytical_velocities': vel_analytical_interp,
                'analytical_true_anomaly': true_anomaly_analytical_interp,
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
        
        # Store results
        self.results[f'single_method_{method_name.lower()}'] = method_results
        self.results['analytical_reference'] = {
            'times': time_ref,
            'positions': pos_analytical,
            'velocities': vel_analytical,
            'true_anomaly': true_anomaly_analytical
        }
        
        return method_results

    def analyze_single_method_results(self, method_name):
        """Analyze and print results for a single method."""
        print(f"\n{'='*80}")
        print(f"{method_name.upper()} METHOD ANALYSIS RESULTS")
        print(f"{'='*80}")
        
        method_key = f'single_method_{method_name.lower()}'
        method_data = self.results[method_key]
        analytical_data = self.results['analytical_reference']
        
        print(f"\nAnalytical Reference Solution:")
        print(f"  Time points: {len(analytical_data['times']):,}")
        
        print(f"\n{method_name} Method Results:")
        print(f"{'='*40}")
        
        # Print summary table
        print(f"{'Step Size':>10} {'Steps':>8} {'Time':>8} {'Max Pos Err':>12} {'Max Vel Err':>12} {'Energy Err':>12} {'True Anom Err':>14}")
        print(f"{'(s)':>10} {'':>8} {'(s)':>8} {'(km)':>12} {'(km/s)':>12} {'(km²/s²)':>12} {'(rad)':>14}")
        print("-" * 94)
        
        for step_size in INTEGRATION_STEPS:
            result = method_data[step_size]
            print(f"{step_size:>10.1f} {result['num_steps']:>8d} "
                  f"{result['computation_time']:>8.3f} "
                  f"{result['max_position_error']:>12.3e} "
                  f"{result['max_velocity_error']:>12.3e} "
                  f"{result['max_energy_error']:>12.3e} "
                  f"{result['max_true_anomaly_error']:>14.3e}")
        
        # Analyze convergence order
        step_sizes = np.array(INTEGRATION_STEPS)
        max_pos_errors = np.array([method_data[h]['max_position_error'] for h in INTEGRATION_STEPS])
        
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

    def save_single_method_results(self, method_name):
        """Save results for a single method."""
        print(f"\nSaving {method_name} results to {RESULTS_DIR}/")
        
        method_key = f'single_method_{method_name.lower()}'
        method_data = self.results[method_key]
        
        # Save summary
        with open(f"{RESULTS_DIR}/{method_name.lower()}_method_summary.txt", 'w') as f:
            f.write(f"Orbit B {method_name} Method Study Results\n")
            f.write("=" * 40 + "\n\n")
            f.write("Orbital Parameters:\n")
            for key, value in self.orbit_params.items():
                f.write(f"  {key}: {value}\n")
            f.write(f"\nIntegration Steps: {INTEGRATION_STEPS}\n")
            f.write(f"Number of Orbits: {NUM_ORBITS_STUDY}\n\n")
            
            f.write(f"{method_name} Method Results:\n")
            f.write("=" * 30 + "\n")
            f.write(f"{'Step Size':>10} {'Steps':>8} {'Time':>8} {'Max Pos Err':>12} {'Max Vel Err':>12} {'True Anom Err':>14}\n")
            f.write(f"{'(s)':>10} {'':>8} {'(s)':>8} {'(km)':>12} {'(km/s)':>12} {'(rad)':>14}\n")
            f.write("-" * 74 + "\n")
            
            for step_size in INTEGRATION_STEPS:
                result = method_data[step_size]
                f.write(f"{step_size:>10.1f} {result['num_steps']:>8d} "
                       f"{result['computation_time']:>8.3f} "
                       f"{result['max_position_error']:>12.3e} "
                       f"{result['max_velocity_error']:>12.3e} "
                       f"{result['max_true_anomaly_error']:>14.3e}\n")
        
        # Save detailed results for each step size
        for step_size in INTEGRATION_STEPS:
            result = method_data[step_size]
            filename = f"{RESULTS_DIR}/{method_name.lower()}_single_method_step_{step_size:g}s.npz"
            
            np.savez(filename,
                     times=result['times'],
                     positions=result['positions'],
                     velocities=result['velocities'],
                     analytical_positions=result['analytical_positions'],
                     analytical_velocities=result['analytical_velocities'],
                     analytical_true_anomaly=result['analytical_true_anomaly'],
                     position_errors=result['position_errors'],
                     velocity_errors=result['velocity_errors'],
                     energy_errors=result['energy_errors'],
                     true_anomaly_method=result['true_anomaly_method'],
                     true_anomaly_errors=result['true_anomaly_errors'],
                     step_size=step_size,
                     method=method_name)
        
        print(f"{method_name} results saved successfully!")

# =============================================================================
# MAIN EXECUTION
# =============================================================================

def run_special_case_only():
    """Run only the special case: RK8 with 100s step for 10 orbits."""
    print("Running SPECIAL CASE ONLY: RK8 with 100s step for 10 orbits...")
    print("=" * 60)
    
    # Create simulation instance
    sim = OrbitSimulation()
    
    # Set up initial conditions
    sim.setup_initial_conditions()
    
    # Run only the special case
    sim.run_special_case_100_orbits()
    
    # Print analysis for special case only
    if 'special_case_100_orbits' in sim.results:
        special_data = sim.results['special_case_100_orbits']
        print(f"\n{'='*60}")
        print("SPECIAL CASE RESULTS SUMMARY")
        print(f"{'='*60}")
        print(f"Method: RK8")
        print(f"Step size: {special_data['step_size']:.1f} s")
        print(f"Number of orbits: {special_data['num_orbits']}")
        print(f"Total steps: {special_data['num_steps']:,}")
        print(f"Computation time: {special_data['computation_time']:.3f} s")
        print(f"Max position error: {special_data['max_position_error']:>12.3e} km")
        print(f"Final position error: {special_data['final_position_error']:>12.3e} km")
        print(f"Max velocity error: {special_data['max_velocity_error']:>12.3e} km/s")
        print(f"Max energy error: {special_data['max_energy_error']:>12.3e} km²/s²")
        
        # Save special case results only
        print(f"\nSaving special case results to {RESULTS_DIR}/")
        Path(RESULTS_DIR).mkdir(exist_ok=True)
        
        filename = f"{RESULTS_DIR}/rk8_special_case_100orbits_100s_standalone.npz"
        np.savez(filename,
                 times=special_data['times'],
                 positions=special_data['positions'],
                 velocities=special_data['velocities'],
                 analytical_positions=special_data['analytical_positions'],
                 analytical_velocities=special_data['analytical_velocities'],
                 position_errors=special_data['position_errors'],
                 velocity_errors=special_data['velocity_errors'],
                 energy_errors=special_data['energy_errors'],
                 step_size=special_data['step_size'],
                 num_orbits=special_data['num_orbits'],
                 method='RK8_special_100_orbits')
        
        # Save summary text file
        with open(f"{RESULTS_DIR}/rk8_special_case_summary.txt", 'w') as f:
            f.write("RK8 Special Case Results: 10 Orbits with 100s Step\n")
            f.write("=" * 50 + "\n\n")
            f.write("Orbital Parameters:\n")
            for key, value in sim.orbit_params.items():
                f.write(f"  {key}: {value}\n")
            f.write(f"\nSimulation Parameters:\n")
            f.write(f"  Method: RK8\n")
            f.write(f"  Step size: {special_data['step_size']:.1f} s\n")
            f.write(f"  Number of orbits: {special_data['num_orbits']}\n")
            f.write(f"  Total steps: {special_data['num_steps']:,}\n")
            f.write(f"\nResults:\n")
            f.write(f"  Computation time: {special_data['computation_time']:.3f} s\n")
            f.write(f"  Max position error: {special_data['max_position_error']:>12.3e} km\n")
            f.write(f"  Final position error: {special_data['final_position_error']:>12.3e} km\n")
            f.write(f"  Max velocity error: {special_data['max_velocity_error']:>12.3e} km/s\n")
            f.write(f"  Max energy error: {special_data['max_energy_error']:>12.3e} km²/s²\n")
        
        print(f"Special case results saved successfully!")
        print(f"  - Numerical data: {filename}")
        print(f"  - Summary: {RESULTS_DIR}/rk8_special_case_summary.txt")
        
        # Create and display plots
        print(f"\nCreating comparison plots...")
        create_special_case_plots(sim)

        # Also create improved plots that explain the eccentric orbit behavior
        print(f"\nCreating IMPROVED plots that explain eccentric orbit physics...")
        create_special_case_plots_improved(sim)
    
    return sim

def run_single_method_study(method_name='RK4'):
    """Run convergence study for a single method (RK4 or RK8) with comprehensive plotting."""
    print(f"Running {method_name} Method Study...")
    print("=" * 60)
    
    # Create simulation instance
    sim = OrbitSimulation()
    
    # Set up initial conditions
    sim.setup_initial_conditions()
    
    # Run single method convergence study
    sim.run_single_method_convergence_study(method_name)
    
    # Analyze results for single method
    sim.analyze_single_method_results(method_name)
    
    # Save results for single method
    sim.save_single_method_results(method_name)
    
    # Create comprehensive plots for single method
    create_single_method_plots(sim, method_name)
    
    print(f"\n{method_name} simulation completed successfully!")
    print(f"Results saved in '{RESULTS_DIR}/' directory")
    
    return sim

def show_menu():
    """Display the interactive menu and get user choice."""
    print("\n" + "="*60)
    print("     ORBIT B NUMERICAL INTEGRATION STUDY")
    print("="*60)
    print()
    print("Select an option:")
    print()
    print("1. Run Full Convergence Study (Both Methods)")
    print("   • Tests RK4 and RK8 with multiple step sizes (1s to 300s)")
    print("   • Analyzes 5 orbits for convergence analysis")
    print("   • Includes special case (RK8, 100s, 10 orbits)")
    print("   • Generates comprehensive plots and analysis")
    print("   • Duration: ~2-5 minutes")
    print()
    print("2. Run RK4 Method Only")
    print("   • Tests RK4 with multiple step sizes (1s to 300s)")
    print("   • Position, velocity, and true anomaly plots vs time")
    print("   • Error analysis plots for each step size")
    print("   • Analyzes 5 orbits for convergence analysis")
    print("   • Duration: ~1-2 minutes")
    print()
    print("3. Run RK8 Method Only")
    print("   • Tests RK8 with multiple step sizes (1s to 300s)")
    print("   • Position, velocity, and true anomaly plots vs time")
    print("   • Error analysis plots for each step size")
    print("   • Analyzes 5 orbits for convergence analysis")
    print("   • Duration: ~1-2 minutes")
    print()
    print("4. Run Special Case Only")
    print("   • RK8 with 100s step size for 10 orbits")
    print("   • Medium-term accuracy assessment (~21.3 days)")
    print("   • Creates comparison plots (position, velocity, true anomaly)")
    print("   • Faster execution, focused results")
    print("   • Duration: ~5-10 seconds")
    print()
    print("5. Help & Information")
    print("   • Show detailed usage instructions")
    print("   • Display orbital parameters")
    print()
    print("6. Exit")
    print()
    print("="*60)
    
    while True:
        try:
            choice = input("Enter your choice (1-6): ").strip()
            if choice in ['1', '2', '3', '4', '5', '6']:
                return int(choice)
            else:
                print("Invalid choice. Please enter 1, 2, 3, 4, 5, or 6.")
        except (KeyboardInterrupt, EOFError):
            print("\nExiting...")
            return 6

def show_help_info():
    """Display detailed help and orbital parameters."""
    print("\n" + "="*60)
    print("     DETAILED INFORMATION")
    print("="*60)
    
    print("\nORBITAL PARAMETERS (Orbit B):")
    print("  • Semi-major axis: 70,000 km")
    print("  • Eccentricity: 0.9 (highly eccentric)")
    print("  • Inclination: 28.5°")
    print("  • Orbital period: ~51.2 hours")
    print("  • Periapsis altitude: ~300 km")
    print("  • Apoapsis altitude: ~133,300 km")
    
    print("\nINTEGRATION METHODS:")
    print("  • RK4: 4th-order Runge-Kutta (classical)")
    print("  • RK8: 8th-order Runge-Kutta (high precision)")
    
    print("\nCONVERGENCE STUDY DETAILS:")
    print("  • Step sizes tested: 1, 5, 10, 30, 60, 100, 300 seconds")
    print("  • Duration: 5 complete orbits (~10.6 days)")
    print("  • Error analysis: Position, velocity, energy, angular momentum")
    print("  • Compares against analytical Kepler solution")
    
    print("\nSPECIAL CASE DETAILS:")
    print("  • Method: RK8 only")
    print("  • Step size: 100 seconds")
    print("  • Duration: 10 complete orbits (~21.3 days)")
    print("  • Purpose: Medium-term accuracy assessment")
    print("  • Direct time-point comparison (no interpolation)")
    
    print("\nOUTPUT FILES:")
    print("  • results/ directory with .npz data files")
    print("  • Convergence summary text files")
    print("  • Error analysis plots (if full study)")
    print("  • Special case comparison plots (position, velocity, true anomaly vs time)")
    
    print("\nCOMMAND LINE OPTIONS:")
    print("  python main_simulation.py                # Interactive menu")
    print("  python main_simulation.py --full         # Skip menu, run full study (both methods)")
    print("  python main_simulation.py --rk4          # Skip menu, run RK4 method only")
    print("  python main_simulation.py --rk8          # Skip menu, run RK8 method only")
    print("  python main_simulation.py --special-case # Skip menu, run special case")
    print("  python main_simulation.py --help         # Show basic help")
    
    print("="*60)
    input("\nPress Enter to return to main menu...")

def run_full_study():
    """Run the complete convergence study."""
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
    
    # Create special case plots if available
    if 'special_case_100_orbits' in sim.results:
        print(f"\nCreating special case comparison plots...")
        create_special_case_plots(sim)
        
        # Also create improved plots that explain the eccentric orbit behavior
        print(f"\nCreating IMPROVED plots that explain eccentric orbit physics...")
        create_special_case_plots_improved(sim)
    
    print("\nSimulation completed successfully!")
    print(f"Results saved in '{RESULTS_DIR}/' directory")

def main():
    """Main function with interactive menu or command line options."""
    import sys
    
    # Handle command line arguments
    if len(sys.argv) > 1:
        arg = sys.argv[1].lower()
        if arg in ["-h", "--help"]:
            show_help_info()
            return
        elif arg == "--special-case":
            run_special_case_only()
            return
        elif arg == "--full":
            run_full_study()
            return
        elif arg == "--rk4":
            run_single_method_study('RK4')
            return
        elif arg == "--rk8":
            run_single_method_study('RK8')
            return
        else:
            print(f"Unknown argument: {sys.argv[1]}")
            print("Use --help to see available options.")
            return
    
    # Interactive menu mode
    while True:
        choice = show_menu()
        
        if choice == 1:
            print("\nStarting Full Convergence Study...")
            run_full_study()
            print(f"\nPress Enter to return to menu...")
            input()
            
        elif choice == 2:
            print("\nStarting RK4 Method Study...")
            run_single_method_study('RK4')
            print(f"\nPress Enter to return to menu...")
            input()
            
        elif choice == 3:
            print("\nStarting RK8 Method Study...")
            run_single_method_study('RK8')
            print(f"\nPress Enter to return to menu...")
            input()
            
        elif choice == 4:
            print("\nStarting Special Case...")
            run_special_case_only()
            print(f"\nPress Enter to return to menu...")
            input()
            
        elif choice == 5:
            show_help_info()
            
        elif choice == 6:
            print("Goodbye!")
            break

def print_usage():
    """Print basic usage instructions."""
    print("Orbit B Numerical Integration Study")
    print("="*40)
    print("Usage:")
    print("  python main_simulation.py           # Interactive menu")
    print("  python main_simulation.py --full    # Run full convergence study (both methods)")
    print("  python main_simulation.py --rk4     # Run RK4 method only")
    print("  python main_simulation.py --rk8     # Run RK8 method only")
    print("  python main_simulation.py --special-case    # Run only the special case")
    print("  python main_simulation.py --help    # Show detailed help")
    print("")
    print("Interactive menu provides detailed options and information.")

if __name__ == "__main__":
    main()