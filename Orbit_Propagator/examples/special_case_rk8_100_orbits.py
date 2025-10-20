"""
Special Case Study: RK8 Integration for 100 Orbits
=================================================

This script runs only the special case simulation: RK8 method with 100s step size
over 100 orbits, without running the full convergence study.

This allows for quick testing and analysis of long-term orbital integration
accuracy without the overhead of the complete convergence analysis.

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
    RESULTS_DIR,
    print_orbit_summary
)
from core.Kepler import (
    analytical_propagation, 
    position_velocity_from_elements
)
from core.RK4 import calculate_orbital_energy, calculate_angular_momentum
from core.RK8 import propagate_orbit_rk8

# =============================================================================
# SPECIAL CASE SIMULATION CLASS
# =============================================================================

class SpecialCaseSimulation:
    """Class for running the RK8 100-orbit special case simulation."""
    
    def __init__(self):
        """Initialize the simulation with Orbit B parameters."""
        self.orbit_params = get_orbit_b_parameters()
        self.results = {}
        
        # Create results directory
        Path(RESULTS_DIR).mkdir(exist_ok=True)
        
        print("RK8 Special Case: 100 Orbits with 100s Step Size")
        print("=" * 55)
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
    
    def run_special_case(self):
        """Run the RK8 special case: 100 orbits with 100s step size."""
        print(f"\n{'='*60}")
        print("RUNNING SPECIAL CASE SIMULATION")
        print(f"{'='*60}")
        print(f"Method: RK8 (8th-order Runge-Kutta)")
        print(f"Step size: 100 seconds")
        print(f"Duration: 100 orbits")
        
        # Parameters for special case
        num_orbits = 100
        step_size = 100.0  # seconds
        period = self.orbit_params['orbital_period']
        t_end = num_orbits * period
        
        print(f"\nSimulation Parameters:")
        print(f"  Orbital period: {period:.1f} s ({period/3600:.2f} hours)")
        print(f"  Total simulation time: {t_end:.1f} s ({t_end/3600:.2f} hours)")
        print(f"  Expected number of steps: ~{int(t_end/step_size):,}")
        
        # Run RK8 integration
        print(f"\nStarting RK8 integration...")
        start_time = time.time()
        times, positions, velocities = propagate_orbit_rk8(
            self.initial_state,
            (0, t_end),
            step_size
        )
        computation_time_rk8 = time.time() - start_time
        
        print(f"  RK8 integration completed!")
        print(f"  Actual steps completed: {len(times):,}")
        print(f"  RK8 computation time: {computation_time_rk8:.3f} s")
        
        # Calculate analytical solution for comparison
        print(f"\nCalculating analytical reference solution...")
        analytical_start = time.time()
        pos_analytical, vel_analytical = analytical_propagation(
            self.orbit_params, times
        )
        computation_time_analytical = time.time() - analytical_start
        print(f"  Analytical computation time: {computation_time_analytical:.3f} s")
        
        # Calculate errors
        print(f"\nCalculating errors and conservation properties...")
        position_errors = np.linalg.norm(positions - pos_analytical, axis=1)
        velocity_errors = np.linalg.norm(velocities - vel_analytical, axis=1)
        
        # Calculate energy conservation
        energy = calculate_orbital_energy(positions, velocities)
        energy_errors = np.abs(energy - energy[0])
        
        # Calculate angular momentum conservation
        angular_momentum = calculate_angular_momentum(positions, velocities)
        momentum_errors = np.abs(angular_momentum - angular_momentum[0])
        
        # Store results
        self.results = {
            'times': times,
            'positions': positions,
            'velocities': velocities,
            'analytical_positions': pos_analytical,
            'analytical_velocities': vel_analytical,
            'position_errors': position_errors,
            'velocity_errors': velocity_errors,
            'energy': energy,
            'energy_errors': energy_errors,
            'angular_momentum': angular_momentum,
            'angular_momentum_errors': momentum_errors,
            'max_position_error': np.max(position_errors),
            'max_velocity_error': np.max(velocity_errors),
            'max_energy_error': np.max(energy_errors),
            'max_momentum_error': np.max(momentum_errors),
            'final_position_error': position_errors[-1],
            'final_velocity_error': velocity_errors[-1],
            'final_energy_error': energy_errors[-1],
            'final_momentum_error': momentum_errors[-1],
            'rk8_computation_time': computation_time_rk8,
            'analytical_computation_time': computation_time_analytical,
            'total_computation_time': computation_time_rk8 + computation_time_analytical,
            'num_steps': len(times),
            'step_size': step_size,
            'num_orbits': num_orbits,
            'simulation_time_span': t_end
        }
        
        return self.results
    
    def analyze_results(self):
        """Analyze and print the simulation results."""
        results = self.results
        
        print(f"\n{'='*70}")
        print("SIMULATION RESULTS ANALYSIS")
        print(f"{'='*70}")
        
        # Basic simulation info
        print(f"\nSimulation Summary:")
        print(f"  Method: RK8 (8th-order Runge-Kutta)")
        print(f"  Step size: {results['step_size']} seconds")
        print(f"  Number of orbits: {results['num_orbits']}")
        print(f"  Total steps: {results['num_steps']:,}")
        print(f"  Simulation time span: {results['simulation_time_span']:.1f} s")
        print(f"  Simulation time span: {results['simulation_time_span']/3600:.2f} hours")
        print(f"  Simulation time span: {results['simulation_time_span']/(3600*24):.2f} days")
        
        # Computation times
        print(f"\nComputation Performance:")
        print(f"  RK8 computation time: {results['rk8_computation_time']:>10.3f} s")
        print(f"  Analytical comp. time: {results['analytical_computation_time']:>10.3f} s")
        print(f"  Total computation time: {results['total_computation_time']:>10.3f} s")
        print(f"  Steps per second (RK8): {results['num_steps']/results['rk8_computation_time']:>10.1f}")
        
        # Error analysis
        print(f"\nAccuracy Analysis:")
        print(f"  Maximum position error: {results['max_position_error']:>12.3e} km")
        print(f"  Final position error: {results['final_position_error']:>12.3e} km")
        print(f"  Maximum velocity error: {results['max_velocity_error']:>12.3e} km/s")
        print(f"  Final velocity error: {results['final_velocity_error']:>12.3e} km/s")
        
        # Conservation analysis
        print(f"\nConservation Analysis:")
        print(f"  Maximum energy error: {results['max_energy_error']:>12.3e} km²/s²")
        print(f"  Final energy error: {results['final_energy_error']:>12.3e} km²/s²")
        print(f"  Max angular momentum error: {results['max_momentum_error']:>12.3e} km²/s")
        print(f"  Final angular momentum error: {results['final_momentum_error']:>12.3e} km²/s")
        
        # Relative error analysis
        initial_energy = results['energy'][0]
        initial_momentum = results['angular_momentum'][0]
        print(f"\nRelative Conservation Errors:")
        print(f"  Relative energy error: {results['final_energy_error']/abs(initial_energy):>12.3e}")
        print(f"  Relative momentum error: {results['final_momentum_error']/abs(initial_momentum):>12.3e}")
        
        # Error growth analysis
        orbit_period = self.orbit_params['orbital_period']
        times_in_orbits = results['times'] / orbit_period
        
        # Find error at different orbital milestones
        milestones = [10, 25, 50, 75, 100]
        print(f"\nError Growth Analysis (Position Error):")
        print(f"{'Orbit':>8} {'Time (h)':>10} {'Position Error (km)':>20}")
        print("-" * 40)
        
        for milestone in milestones:
            if milestone <= results['num_orbits']:
                # Find closest index to this orbital milestone
                target_time = milestone * orbit_period
                idx = np.argmin(np.abs(results['times'] - target_time))
                actual_orbit = results['times'][idx] / orbit_period
                actual_time_hours = results['times'][idx] / 3600
                pos_error = results['position_errors'][idx]
                print(f"{actual_orbit:>8.1f} {actual_time_hours:>10.1f} {pos_error:>20.3e}")
    
    def save_results(self):
        """Save results to files."""
        print(f"\nSaving results to {RESULTS_DIR}/")
        
        results = self.results
        
        # Save detailed numerical results
        filename = f"{RESULTS_DIR}/rk8_special_case_100orbits_100s_standalone.npz"
        np.savez(filename,
                 times=results['times'],
                 positions=results['positions'],
                 velocities=results['velocities'],
                 analytical_positions=results['analytical_positions'],
                 analytical_velocities=results['analytical_velocities'],
                 position_errors=results['position_errors'],
                 velocity_errors=results['velocity_errors'],
                 energy=results['energy'],
                 energy_errors=results['energy_errors'],
                 angular_momentum=results['angular_momentum'],
                 angular_momentum_errors=results['angular_momentum_errors'],
                 step_size=results['step_size'],
                 num_orbits=results['num_orbits'],
                 method='RK8_special_100_orbits_standalone')
        
        # Save summary report
        summary_filename = f"{RESULTS_DIR}/rk8_special_case_summary.txt"
        with open(summary_filename, 'w') as f:
            f.write("RK8 Special Case: 100 Orbits Simulation Summary\n")
            f.write("=" * 50 + "\n\n")
            
            f.write("Orbital Parameters:\n")
            for key, value in self.orbit_params.items():
                f.write(f"  {key}: {value}\n")
            
            f.write(f"\nSimulation Parameters:\n")
            f.write(f"  Method: RK8 (8th-order Runge-Kutta)\n")
            f.write(f"  Step size: {results['step_size']} seconds\n")
            f.write(f"  Number of orbits: {results['num_orbits']}\n")
            f.write(f"  Total steps: {results['num_steps']:,}\n")
            f.write(f"  Simulation time span: {results['simulation_time_span']:.1f} s\n")
            
            f.write(f"\nPerformance Results:\n")
            f.write(f"  RK8 computation time: {results['rk8_computation_time']:.3f} s\n")
            f.write(f"  Analytical computation time: {results['analytical_computation_time']:.3f} s\n")
            f.write(f"  Total computation time: {results['total_computation_time']:.3f} s\n")
            
            f.write(f"\nAccuracy Results:\n")
            f.write(f"  Max position error: {results['max_position_error']:.3e} km\n")
            f.write(f"  Final position error: {results['final_position_error']:.3e} km\n")
            f.write(f"  Max velocity error: {results['max_velocity_error']:.3e} km/s\n")
            f.write(f"  Final velocity error: {results['final_velocity_error']:.3e} km/s\n")
            f.write(f"  Max energy error: {results['max_energy_error']:.3e} km²/s²\n")
            f.write(f"  Final energy error: {results['final_energy_error']:.3e} km²/s²\n")
        
        print(f"  Numerical data saved to: {filename}")
        print(f"  Summary report saved to: {summary_filename}")
        print(f"Results saved successfully!")
    
    def create_position_error_plot(self):
        """Create and save the position error plot."""
        print(f"\nCreating position error plot...")
        
        results = self.results
        
        fig = plt.figure(figsize=(16, 10))
        ax = fig.add_subplot(1, 1, 1)
        
        # Convert time to orbital periods
        orbit_period = self.orbit_params['orbital_period']
        times_orbits = results['times'] / orbit_period
        position_errors = results['position_errors']
        
        # Plot position error
        ax.semilogy(times_orbits, position_errors, 
                   'b-', linewidth=2, alpha=0.8, label='RK8 Position Error (Δt=100s)')
        
        # Add reference lines
        max_error = results['max_position_error']
        final_error = results['final_position_error']
        mean_error = np.mean(position_errors)
        
        ax.axhline(y=max_error, color='red', linestyle='--', alpha=0.7, 
                   label=f'Maximum Error: {max_error:.2e} km')
        ax.axhline(y=final_error, color='orange', linestyle='--', alpha=0.7, 
                   label=f'Final Error: {final_error:.2e} km')
        ax.axhline(y=mean_error, color='green', linestyle='--', alpha=0.7, 
                   label=f'Mean Error: {mean_error:.2e} km')
        
        # Add vertical grid lines every 10 orbits
        for orbit_mark in range(10, 101, 10):
            ax.axvline(x=orbit_mark, color='gray', linestyle=':', alpha=0.3)
            if orbit_mark % 20 == 0:
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
        ax.set_xlim(0, 100)
        
        # Add info box
        info_text = (f"Integration Method: RK8\n"
                    f"Step Size: {results['step_size']} s\n"
                    f"Total Steps: {results['num_steps']:,}\n"
                    f"RK8 Computation Time: {results['rk8_computation_time']:.1f} s\n"
                    f"Simulation Span: {results['num_orbits']} orbits")
        
        ax.text(0.02, 0.98, info_text, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        # Save plot
        plot_filename = f'{RESULTS_DIR}/rk8_special_case_100orbits_position_error_standalone.png'
        plt.tight_layout()
        plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"  Position error plot saved to: {plot_filename}")

# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main():
    """Run the RK8 special case simulation."""
    print("Starting RK8 Special Case: 100 Orbits with 100s Step Size")
    print(f"Current working directory: {os.getcwd()}")
    
    # Create simulation object
    sim = SpecialCaseSimulation()
    
    # Set up initial conditions
    sim.setup_initial_conditions()
    
    # Run the special case
    sim.run_special_case()
    
    # Analyze results
    sim.analyze_results()
    
    # Save results
    sim.save_results()
    
    # Create plot
    sim.create_position_error_plot()
    
    print(f"\n{'='*70}")
    print("SPECIAL CASE SIMULATION COMPLETED SUCCESSFULLY!")
    print(f"{'='*70}")
    print(f"Results saved in '{RESULTS_DIR}/' directory")
    print(f"Files created:")
    print(f"  - rk8_special_case_100orbits_100s_standalone.npz")
    print(f"  - rk8_special_case_summary.txt")
    print(f"  - rk8_special_case_100orbits_position_error_standalone.png")

if __name__ == "__main__":
    main()