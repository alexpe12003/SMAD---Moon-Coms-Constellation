"""
Porkchop Plot Analysis with 150-Hour Simulation
==============================================
Comprehensive porkchop plot generation using 150 hours of simulated orbit data.

Grid Parameters:
- X-axis: Departure time (0 to 2 LEO satellite orbits)
- Y-axis: Time of flight (10 to 140 hours)
- Uses interpolated position data from 150-hour simulation
- Lambert solver calculates transfer orbits and delta-V requirements
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import math
import sys
import os

# Add the lambert solver to the path (go up one more level since we're in core_modules/)
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'Duda_Lambert', 'lambert')))

import lambert
from porkchop import create_leo_satellite, create_moon_orbit

# Physical Constants
MU_EARTH_KM3_S2 = 398600.4418

def quick_extended_porkchop():
    """Generate porkchop plot using 120 hours of simulated orbit data"""
    
    print("PORKCHOP ANALYSIS WITH 150-HOUR SIMULATION DATA")
    print("=" * 50)
    
    # Create orbit objects
    leo_satellite = create_leo_satellite()
    moon = create_moon_orbit()
    
    print(f"LEO period: {leo_satellite.period/60:.1f} minutes")
    print(f"Moon period: {moon.period/3600/24:.1f} days")
    
    # Define parameters using 150-hour simulation data
    dep_duration = 2 * leo_satellite.period    # 2 LEO orbits (x-axis)
    tof_min = 10 * 3600                        # 10 hours minimum (y-axis)
    tof_max = 140 * 3600                       # 140 hours maximum (y-axis)
    
    # Grid resolution - balanced for good resolution and computation time
    n_dep = 500    # 500 departure points across 2 orbits
    n_tof = 500    # 500 time of flight points from 10-140 hours
    
    departure_times = np.linspace(0, dep_duration, n_dep)
    tof_times = np.linspace(tof_min, tof_max, n_tof)
    
    print(f"Analysis grid: {n_dep} x {n_tof} = {n_dep * n_tof} calculations")
    print(f"Departure window: 0 to {dep_duration/3600:.1f} hours ({dep_duration/leo_satellite.period:.1f} LEO orbits)")
    print(f"TOF range: {tof_min/3600:.0f} to {tof_max/3600:.0f} hours")
    
    # Verify we stay within 150-hour simulation bounds
    max_mission_time = dep_duration + tof_max
    print(f"Maximum mission time: {max_mission_time/3600:.1f} hours")
    if max_mission_time <= 150 * 3600:
        print("✅ All scenarios within 150-hour simulation bounds")
    else:
        print("⚠️  Some scenarios exceed 150-hour simulation bounds")
    
    # Generate 150 hours of orbit data with 1-minute intervals
    print("Generating 150 hours of orbit data...")
    simulation_hours = 150
    leo_orbits_needed = simulation_hours / (leo_satellite.period / 3600)
    moon_orbits_needed = simulation_hours / (moon.period / 3600)
    
    leo_data = leo_satellite.generate_orbit_data(time_interval_s=60, num_orbits=leo_orbits_needed)
    moon_data = moon.generate_orbit_data(time_interval_s=60, num_orbits=moon_orbits_needed)
    
    print(f"Generated {len(leo_data['time_s'])} data points for each orbit")
    
    def get_position_from_data(data, time_s):
        """Interpolate position from orbit data at given time"""
        # Find closest data point (linear interpolation for better accuracy)
        time_array = data['time_s']
        if time_s <= time_array[0]:
            idx = 0
        elif time_s >= time_array[-1]:
            idx = -1
        else:
            idx = np.searchsorted(time_array, time_s)
            if idx > 0:
                # Linear interpolation between two points
                t1, t2 = time_array[idx-1], time_array[idx]
                alpha = (time_s - t1) / (t2 - t1)
                x = data['x_km'][idx-1] * (1 - alpha) + data['x_km'][idx] * alpha
                y = data['y_km'][idx-1] * (1 - alpha) + data['y_km'][idx] * alpha
                # Calculate velocity using finite difference
                vx = (data['x_km'][idx] - data['x_km'][idx-1]) / (t2 - t1)
                vy = (data['y_km'][idx] - data['y_km'][idx-1]) / (t2 - t1)
                return np.array([x, y, 0]), np.array([vx, vy, 0])
        
        # Use exact data point
        x, y = data['x_km'][idx], data['y_km'][idx]
        # Calculate velocity from orbital mechanics
        radius = np.sqrt(x**2 + y**2)
        if 'leo' in str(type(data)).lower() or radius < 50000:  # LEO satellite
            velocity_mag = leo_satellite.velocity
        else:  # Moon
            velocity_mag = moon.velocity
        
        # Velocity is perpendicular to position vector
        vx = -velocity_mag * y / radius
        vy = velocity_mag * x / radius
        
        return np.array([x, y, 0]), np.array([vx, vy, 0])
    
    # Initialize results
    dv_total = np.full((n_tof, n_dep), np.nan)
    dv_departure = np.full((n_tof, n_dep), np.nan)
    dv_arrival = np.full((n_tof, n_dep), np.nan)
    
    print("Calculating transfer solutions using 120-hour simulation data...")
    valid_count = 0
    
    for i, dep_time in enumerate(departure_times):
        progress = (i + 1) / n_dep * 100
        print(f"Progress: {progress:.1f}% (Column {i+1}/{n_dep})")
        
        for j, tof in enumerate(tof_times):
            try:
                # Check if arrival time is within simulation bounds
                arrival_time = dep_time + tof
                if arrival_time > 150 * 3600:  # Beyond 150 hours
                    continue
                
                # Get satellite state at departure from simulation data
                R1, V1_circ = get_position_from_data(leo_data, dep_time)
                
                # Get Moon state at arrival from simulation data
                R2, V2_circ = get_position_from_data(moon_data, arrival_time)
                
                # Solve Lambert problem
                V1_trans, V2_trans, _, _ = lambert.lambert(R1, R2, tof, 'pro', MU_EARTH_KM3_S2)
                
                # Calculate delta-V
                dv_dep = np.linalg.norm(V1_trans - V1_circ)
                dv_arr = np.linalg.norm(V2_trans - V2_circ)
                dv_tot = dv_dep + dv_arr
                
                # Store if reasonable (filter extreme values)
                if dv_tot < 20:  # Increased threshold for longer TOF range
                    dv_total[j, i] = dv_tot
                    dv_departure[j, i] = dv_dep
                    dv_arrival[j, i] = dv_arr
                    valid_count += 1
                
            except Exception as e:
                # Skip failed calculations
                continue
    
    print(f"Valid solutions found: {valid_count}/{n_dep * n_tof}")
    
    # Create plots
    dep_grid, tof_grid = np.meshgrid(departure_times / 3600, tof_times / 3600)  # Convert to hours
    
    def save_plot(data, title, filename):
        """Save a porkchop plot"""
        fig, ax = plt.subplots(figsize=(16, 12))
        
        if not np.all(np.isnan(data)):
            # Contour plot
            levels = 25
            contour = ax.contourf(dep_grid, tof_grid, data, levels=levels, cmap='plasma')
            
            # Contour lines
            contour_lines = ax.contour(dep_grid, tof_grid, data, levels=15, 
                                     colors='white', alpha=0.7, linewidths=1)
            ax.clabel(contour_lines, inline=True, fontsize=10, fmt='%.2f')
            
            # Colorbar
            cbar = plt.colorbar(contour, ax=ax)
            cbar.set_label('ΔV (km/s)', rotation=270, labelpad=20, fontsize=12)
            
            # Find optimal
            min_idx = np.unravel_index(np.nanargmin(data), data.shape)
            opt_dep = dep_grid[min_idx]
            opt_tof = tof_grid[min_idx]
            opt_dv = data[min_idx]
            
            # Mark optimal point
            ax.plot(opt_dep, opt_tof, 'r*', markersize=25, 
                   markeredgecolor='white', markeredgewidth=3)
            
            # LEO orbit lines - mark each complete orbit
            leo_period_hours = leo_satellite.period / 3600
            for n in range(1, 3):  # Up to 2 orbits
                x_pos = n * leo_period_hours
                if x_pos <= np.max(dep_grid):
                    ax.axvline(x_pos, color='cyan', linestyle='--', alpha=0.8, linewidth=2)
                    ax.text(x_pos + 0.05, ax.get_ylim()[1] * 0.95, 
                           f'{n} orbit', rotation=90, color='cyan', fontsize=10)
            
            # Add time markers on TOF axis
            for tof_mark in [24, 48, 72, 96, 120]:  # Days markers up to 5 days
                if tof_mark <= np.max(tof_grid):
                    ax.axhline(tof_mark, color='orange', linestyle=':', alpha=0.6)
                    ax.text(ax.get_xlim()[1] * 0.95, tof_mark + 1, 
                           f'{tof_mark/24:.0f}d', color='orange', fontsize=10)
            
            # Stats box
            stats_text = f'Optimal Solution:\\nDeparture: {opt_dep:.2f} h\\nTOF: {opt_tof:.1f} h ({opt_tof/24:.1f} days)\\nΔV: {opt_dv:.3f} km/s'
            ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
                   verticalalignment='top', fontsize=12,
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
            
            print(f"Optimal for {filename}: ΔV={opt_dv:.3f} km/s at {opt_dep:.2f} h departure, {opt_tof:.1f} h TOF")
        
        ax.set_xlabel('Departure Time (hours from simulation start)', fontsize=14)
        ax.set_ylabel('Time of Flight (hours)', fontsize=14)
        ax.set_title(title, fontsize=16, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved: {filename}")
    
    # Generate plots - save to final_results directory
    results_dir = os.path.join(os.path.dirname(__file__), '..', 'final_results')
    print("\\nCreating porkchop plots...")
    save_plot(dv_total, 'Porkchop Plot: Total ΔV (150h Simulation)', os.path.join(results_dir, 'porkchop_150h_total.png'))
    save_plot(dv_departure, 'Porkchop Plot: Departure ΔV (150h Simulation)', os.path.join(results_dir, 'porkchop_150h_departure.png'))
    save_plot(dv_arrival, 'Porkchop Plot: Arrival ΔV (150h Simulation)', os.path.join(results_dir, 'porkchop_150h_arrival.png'))
    
    # Final summary
    if not np.all(np.isnan(dv_total)):
        min_idx = np.unravel_index(np.nanargmin(dv_total), dv_total.shape)
        
        print(f"\\n🎯 OPTIMAL SOLUTION FROM 150-HOUR SIMULATION:")
        print(f"   Departure time: {dep_grid[min_idx]:.2f} hours ({dep_grid[min_idx]/leo_satellite.period*3600:.1f} LEO orbits)")
        print(f"   Time of flight: {tof_grid[min_idx]:.1f} hours ({tof_grid[min_idx]/24:.1f} days)")
        print(f"   Total ΔV: {dv_total[min_idx]:.3f} km/s")
        print(f"   Departure ΔV: {dv_departure[min_idx]:.3f} km/s")
        print(f"   Arrival ΔV: {dv_arrival[min_idx]:.3f} km/s")
        
        # Mission timeline
        dep_time_actual = departure_times[min_idx[1]]
        tof_actual = tof_times[min_idx[0]]
        arrival_time_actual = dep_time_actual + tof_actual
        
        print(f"\\n📊 MISSION TIMELINE:")
        print(f"   Simulation start: T+0 hours")
        print(f"   Departure: T+{dep_time_actual/3600:.2f} hours")
        print(f"   Arrival at Moon: T+{arrival_time_actual/3600:.1f} hours")
        print(f"   Total mission time: {arrival_time_actual/3600:.1f} hours ({arrival_time_actual/3600/24:.1f} days)")
        
        # Simulation verification
        simulation_coverage = arrival_time_actual / (150 * 3600) * 100
        print(f"\\n✅ SIMULATION VERIFICATION:")
        print(f"   Mission uses {simulation_coverage:.1f}% of 150-hour simulation")
        print(f"   Departure window: 0 to {dep_duration/3600:.1f} hours ({dep_duration/leo_satellite.period:.1f} LEO orbits)")
        print(f"   TOF range: {tof_min/3600:.0f} to {tof_max/3600:.0f} hours")
        print(f"   All scenarios within simulation bounds: {'Yes' if max_mission_time <= 150*3600 else 'No'}")
        
        # Performance statistics
        total_calculations = n_dep * n_tof
        success_rate = valid_count / total_calculations * 100
        print(f"\\n📈 ANALYSIS STATISTICS:")
        print(f"   Grid size: {n_dep} × {n_tof} = {total_calculations} calculations")
        print(f"   Valid solutions: {valid_count} ({success_rate:.1f}%)")
        print(f"   Resolution: {dep_duration/n_dep/60:.1f} min departure, {(tof_max-tof_min)/n_tof/3600:.1f} h TOF")
    
    return dv_total, dv_departure, dv_arrival

if __name__ == "__main__":
    results = quick_extended_porkchop()
    print("\\n✅ 150-hour simulation porkchop analysis complete!")