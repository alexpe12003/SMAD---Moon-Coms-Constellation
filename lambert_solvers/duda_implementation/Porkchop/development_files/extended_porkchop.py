"""
Extended Porkchop Plot Generator with Proper Moon Simulation
===========================================================
Generates porkchop plots with extended Moon orbital data to properly cover
all departure times and time of flights.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import math
import sys
import os

# Add the lambert solver to the path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'Duda_Lambert', 'lambert')))

import lambert

# Physical Constants
MU_EARTH_KM3_S2 = 398600.4418
EARTH_RADIUS_KM = 6378
MOON_DISTANCE_KM = 384400

class CircularOrbit:
    """Class to represent a circular orbit and calculate positions over time"""
    
    def __init__(self, radius_km, mu_km3_s2, name="Orbit"):
        self.radius = radius_km
        self.mu = mu_km3_s2
        self.name = name
        
        # Calculate orbital parameters
        self.velocity = math.sqrt(mu_km3_s2 / radius_km)
        self.period = 2 * math.pi * math.sqrt(radius_km**3 / mu_km3_s2)
        self.angular_velocity = 2 * math.pi / self.period
        
        print(f"\n{self.name} Orbital Parameters:")
        print(f"  Radius: {self.radius} km")
        print(f"  Velocity: {self.velocity:.3f} km/s")
        print(f"  Period: {self.period/3600:.2f} hours ({self.period/60:.1f} minutes)")
        print(f"  Angular velocity: {self.angular_velocity:.6f} rad/s")
        print(f"  Angular velocity: {math.degrees(self.angular_velocity):.6f} deg/s")
    
    def get_position_at_time(self, time_s):
        """Calculate position at given time in seconds from initial position"""
        theta_rad = self.angular_velocity * time_s
        theta_deg = math.degrees(theta_rad)
        
        x = self.radius * math.cos(theta_rad)
        y = self.radius * math.sin(theta_rad)
        
        return theta_rad, theta_deg, x, y

def create_leo_satellite():
    """Create LEO satellite orbit (200km altitude)"""
    leo_radius = EARTH_RADIUS_KM + 200  # 200km altitude
    return CircularOrbit(leo_radius, MU_EARTH_KM3_S2, "LEO Satellite (200km)")

def create_moon_orbit():
    """Create Moon orbit (average Earth-Moon distance)"""
    return CircularOrbit(MOON_DISTANCE_KM, MU_EARTH_KM3_S2, "Moon")

def generate_extended_porkchop_plot():
    """Generate porkchop plot with proper extended Moon simulation"""
    
    print("EXTENDED PORKCHOP PLOT WITH PROPER MOON SIMULATION")
    print("=" * 55)
    
    # Create orbit objects
    leo_satellite = create_leo_satellite()
    moon = create_moon_orbit()
    
    print(f"\nLEO period: {leo_satellite.period/60:.1f} minutes")
    print(f"Moon period: {moon.period/3600/24:.1f} days")
    
    # Define parameters for comprehensive analysis
    dep_duration = 3 * leo_satellite.period  # 3 LEO orbits for departure window
    tof_min = 18 * 3600                      # 18 hours minimum TOF
    tof_max = 120 * 3600                     # 5 days maximum TOF
    
    # Calculate total simulation time needed
    total_sim_time = dep_duration + tof_max + 7200  # Add 2h buffer
    total_sim_hours = total_sim_time / 3600
    
    print(f"Analysis parameters:")
    print(f"  Departure window: {dep_duration/3600:.1f} hours ({dep_duration/leo_satellite.period:.1f} LEO orbits)")
    print(f"  Time of flight range: {tof_min/3600:.1f} to {tof_max/3600:.1f} hours")
    print(f"  Total simulation time required: {total_sim_hours:.1f} hours ({total_sim_hours/24:.1f} days)")
    
    # Create analysis grids - using higher resolution
    departure_times = np.linspace(0, dep_duration, 40)  # 40 departure points
    tof_times = np.linspace(tof_min, tof_max, 35)       # 35 time of flight points
    
    print(f"Computational grid: {len(departure_times)} x {len(tof_times)} = {len(departure_times) * len(tof_times)} calculations")
    
    # Pre-calculate Moon positions with finer resolution for accuracy
    moon_time_step = 900  # 15 minute intervals for Moon positions
    moon_times = np.arange(0, total_sim_time + moon_time_step, moon_time_step)
    
    print(f"Pre-calculating {len(moon_times)} Moon positions over {total_sim_hours:.1f} hours...")
    moon_positions = []
    moon_velocities = []
    
    for i, t in enumerate(moon_times):
        if i % 100 == 0:
            progress = i / len(moon_times) * 100
            print(f"  Moon calculation progress: {progress:.1f}%")
        
        _, _, x, y = moon.get_position_at_time(t)
        moon_positions.append([x, y, 0])
        
        # Calculate velocity for Moon at this time (for delta-V calculation)
        theta = moon.angular_velocity * t
        vx = -moon.velocity * math.sin(theta)
        vy = moon.velocity * math.cos(theta)
        moon_velocities.append([vx, vy, 0])
    
    moon_positions = np.array(moon_positions)
    moon_velocities = np.array(moon_velocities)
    print("Moon position calculation complete!")
    
    # Initialize result matrices
    dv_total = np.full((len(tof_times), len(departure_times)), np.nan)
    dv_departure = np.full((len(tof_times), len(departure_times)), np.nan)
    dv_arrival = np.full((len(tof_times), len(departure_times)), np.nan)
    
    print(f"\nCalculating {len(departure_times) * len(tof_times)} transfer solutions...")
    
    # Track calculation statistics
    valid_solutions = 0
    failed_lambert = 0
    high_dv_filtered = 0
    
    for i, dep_time in enumerate(departure_times):
        # Progress update
        progress = (i + 1) / len(departure_times) * 100
        print(f"Row {i+1}/{len(departure_times)} ({progress:.1f}%) - Dep time: {dep_time/60:.1f} min")
        
        for j, tof in enumerate(tof_times):
            try:
                # Get LEO position and velocity at departure
                dep_theta = leo_satellite.angular_velocity * dep_time
                x1 = leo_satellite.radius * math.cos(dep_theta)
                y1 = leo_satellite.radius * math.sin(dep_theta)
                R1 = np.array([x1, y1, 0])
                
                # LEO velocity vector at departure
                v1_x = -leo_satellite.velocity * math.sin(dep_theta)
                v1_y = leo_satellite.velocity * math.cos(dep_theta)
                V1_circ = np.array([v1_x, v1_y, 0])
                
                # Get Moon position at arrival (interpolate from pre-calculated data)
                arrival_time = dep_time + tof
                idx = np.searchsorted(moon_times, arrival_time)
                
                if idx == 0:
                    R2 = moon_positions[0]
                    V2_circ = moon_velocities[0]
                elif idx >= len(moon_positions):
                    R2 = moon_positions[-1]
                    V2_circ = moon_velocities[-1]
                else:
                    # Linear interpolation for position and velocity
                    t1, t2 = moon_times[idx-1], moon_times[idx]
                    alpha = (arrival_time - t1) / (t2 - t1)
                    R2 = moon_positions[idx-1] + alpha * (moon_positions[idx] - moon_positions[idx-1])
                    V2_circ = moon_velocities[idx-1] + alpha * (moon_velocities[idx] - moon_velocities[idx-1])
                
                # Solve Lambert's problem
                V1_trans, V2_trans, _, _ = lambert.lambert(R1, R2, tof, 'pro', MU_EARTH_KM3_S2)
                
                # Calculate delta-V components
                dv_dep = np.linalg.norm(V1_trans - V1_circ)
                dv_arr = np.linalg.norm(V2_trans - V2_circ)
                dv_tot = dv_dep + dv_arr
                
                # Store if reasonable (filter out unrealistic solutions)
                if dv_tot < 20:  # Upper limit for realistic transfers
                    dv_total[j, i] = dv_tot
                    dv_departure[j, i] = dv_dep
                    dv_arrival[j, i] = dv_arr
                    valid_solutions += 1
                else:
                    high_dv_filtered += 1
                
            except Exception as e:
                failed_lambert += 1
                continue
    
    print(f"\nCalculation Statistics:")
    print(f"  Valid solutions: {valid_solutions}")
    print(f"  Failed Lambert solutions: {failed_lambert}")
    print(f"  High ΔV filtered: {high_dv_filtered}")
    print(f"  Total calculations: {len(departure_times) * len(tof_times)}")
    
    # Convert to plotting coordinates
    dep_grid, tof_grid = np.meshgrid(departure_times / 60, tof_times / 3600)
    
    # Create enhanced plots
    def create_enhanced_plot(data_matrix, title, filename, units="km/s"):
        """Create and save enhanced porkchop plot"""
        fig, ax = plt.subplots(figsize=(16, 12))
        
        valid_mask = ~np.isnan(data_matrix)
        
        if np.any(valid_mask):
            # Enhanced color scale
            vmin = np.nanpercentile(data_matrix, 1)
            vmax = np.nanpercentile(data_matrix, 99)
            
            # High-resolution contour plot
            levels = np.linspace(vmin, vmax, 30)
            contour = ax.contourf(dep_grid, tof_grid, data_matrix, 
                                levels=levels, cmap='viridis', extend='both')
            
            # Detailed contour lines with labels
            contour_lines = ax.contour(dep_grid, tof_grid, data_matrix, 
                                     levels=20, colors='white', alpha=0.8, linewidths=0.6)
            ax.clabel(contour_lines, inline=True, fontsize=9, fmt='%.2f')
            
            # Enhanced colorbar
            cbar = plt.colorbar(contour, ax=ax, shrink=0.8)
            cbar.set_label(f'ΔV ({units})', rotation=270, labelpad=25, fontsize=14)
            cbar.ax.tick_params(labelsize=12)
            
            # Find and mark optimal solution
            min_idx = np.unravel_index(np.nanargmin(data_matrix), data_matrix.shape)
            opt_dep = dep_grid[min_idx]
            opt_tof = tof_grid[min_idx]
            opt_dv = data_matrix[min_idx]
            
            ax.plot(opt_dep, opt_tof, 'r*', markersize=25, 
                   markeredgecolor='white', markeredgewidth=3,
                   label=f'Optimal: {opt_dv:.3f} {units}')
            
            # Add LEO orbit reference lines
            leo_period_min = leo_satellite.period / 60
            for n in range(1, 5):
                x_pos = n * leo_period_min
                if x_pos <= np.max(dep_grid):
                    ax.axvline(x_pos, color='cyan', linestyle='--', 
                              alpha=0.9, linewidth=1.5)
                    ax.text(x_pos + 3, ax.get_ylim()[1] * 0.97, 
                           f'{n} LEO orbit', rotation=90, color='cyan', 
                           fontsize=11, fontweight='bold')
            
            # Enhanced statistics box
            leo_orbit_fraction = (opt_dep % leo_period_min) / leo_period_min
            stats_text = (f'🎯 OPTIMAL SOLUTION:\n'
                         f'Departure: {opt_dep:.1f} min\n'
                         f'Time of Flight: {opt_tof:.1f} h ({opt_tof/24:.1f} days)\n'
                         f'ΔV: {opt_dv:.3f} {units}\n'
                         f'LEO phase: {leo_orbit_fraction*100:.1f}%')
            
            ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
                   verticalalignment='top', fontsize=12, fontweight='bold',
                   bbox=dict(boxstyle='round,pad=0.7', facecolor='lightblue', 
                            alpha=0.95, edgecolor='navy'))
            
            print(f"📊 {filename}: Optimal ΔV = {opt_dv:.3f} {units} at {opt_dep:.1f} min dep, {opt_tof:.1f} h TOF")
            
        else:
            ax.text(0.5, 0.5, 'No valid solutions found', 
                   transform=ax.transAxes, ha='center', va='center', 
                   fontsize=18, color='red')
        
        # Enhanced formatting
        ax.set_xlabel('Departure Time (minutes from initial LEO position)', fontsize=14)
        ax.set_ylabel('Time of Flight (hours)', fontsize=14)
        ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
        ax.grid(True, alpha=0.4, linestyle='-', linewidth=0.5)
        if 'opt_dv' in locals():
            ax.legend(fontsize=12, loc='upper right')
        
        # Improved tick formatting
        ax.tick_params(labelsize=12)
        
        plt.tight_layout()
        plt.savefig(filename, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        print(f"✅ Saved: {filename}")
    
    # Generate all plots
    print(f"\n🎨 Creating enhanced porkchop plots...")
    create_enhanced_plot(dv_total, 'LEO-to-Moon Transfer: Total ΔV Requirements', 
                        'extended_porkchop_total_deltav.png')
    create_enhanced_plot(dv_departure, 'LEO-to-Moon Transfer: Departure ΔV Component', 
                        'extended_porkchop_departure_deltav.png')
    create_enhanced_plot(dv_arrival, 'LEO-to-Moon Transfer: Arrival ΔV Component', 
                        'extended_porkchop_arrival_deltav.png')
    
    # Comprehensive final analysis
    if not np.all(np.isnan(dv_total)):
        min_idx = np.unravel_index(np.nanargmin(dv_total), dv_total.shape)
        
        opt_dep_time = dep_grid[min_idx]
        opt_tof = tof_grid[min_idx]
        opt_dv_total = dv_total[min_idx]
        opt_dv_dep = dv_departure[min_idx]
        opt_dv_arr = dv_arrival[min_idx]
        
        print(f"\n🚀 FINAL COMPREHENSIVE ANALYSIS:")
        print(f"{'='*50}")
        print(f"🎯 OPTIMAL TRANSFER SOLUTION:")
        print(f"   Departure time: {opt_dep_time:.1f} minutes ({opt_dep_time/60:.2f} hours)")
        print(f"   Time of flight: {opt_tof:.1f} hours ({opt_tof/24:.2f} days)")
        print(f"   Total ΔV: {opt_dv_total:.3f} km/s")
        print(f"   • Departure ΔV: {opt_dv_dep:.3f} km/s ({opt_dv_dep/opt_dv_total*100:.1f}%)")
        print(f"   • Arrival ΔV: {opt_dv_arr:.3f} km/s ({opt_dv_arr/opt_dv_total*100:.1f}%)")
        
        # Additional orbital mechanics insights
        leo_orbit_num = opt_dep_time / (leo_satellite.period/60)
        leo_orbit_fraction = (opt_dep_time % (leo_satellite.period/60)) / (leo_satellite.period/60)
        
        print(f"\n📐 ORBITAL MECHANICS INSIGHTS:")
        print(f"   LEO orbit number at departure: {leo_orbit_num:.2f}")
        print(f"   Phase in LEO orbit: {leo_orbit_fraction*100:.1f}% ({leo_orbit_fraction*360:.1f}°)")
        print(f"   Energy efficiency: {(1-opt_dv_arr/opt_dv_total)*100:.1f}% front-loaded")
        
        # Mission planning recommendations
        print(f"\n🛰️ MISSION PLANNING RECOMMENDATIONS:")
        if opt_dv_arr < 0.1:
            print(f"   ✅ Excellent arrival efficiency - minimal lunar insertion burn required")
        if opt_dv_dep > 3.0:
            print(f"   ⚠️  High departure ΔV - consider propulsion system requirements")
        if opt_tof < 72:
            print(f"   ⏱️  Fast transfer - good for crew missions")
        else:
            print(f"   🐌 Slow transfer - suitable for cargo missions")
    
    return dv_total, dv_departure, dv_arrival, dep_grid, tof_grid

if __name__ == "__main__":
    print("🌙 Extended Porkchop Plot Analysis Starting...")
    results = generate_extended_porkchop_plot()
    print("\n🎉 Extended porkchop analysis complete! Check the generated PNG files.")