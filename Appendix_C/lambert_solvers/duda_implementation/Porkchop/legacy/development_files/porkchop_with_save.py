"""
File-Saving Porkchop Plot Generator
===================================
Generates porkchop plots and saves them to image files.
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
from porkchop import create_leo_satellite, create_moon_orbit

# Physical Constants
MU_EARTH_KM3_S2 = 398600.4418

def generate_porkchop_plot_with_save():
    """Generate porkchop plot and save to file"""
    
    print("PORKCHOP PLOT WITH FILE SAVE")
    print("=" * 40)
    
    # Create orbit objects
    leo_satellite = create_leo_satellite()
    moon = create_moon_orbit()
    
    print(f"LEO period: {leo_satellite.period/60:.1f} minutes")
    
    # Define parameters for good resolution but reasonable speed
    dep_duration = 3 * leo_satellite.period  # 3 LEO orbits
    tof_min = 18 * 3600                      # 18 hours minimum
    tof_max = 96 * 3600                      # 4 days maximum
    
    # Create grids
    departure_times = np.linspace(0, dep_duration, 30)  # 30 departure points
    tof_times = np.linspace(tof_min, tof_max, 25)       # 25 time of flight points
    
    print(f"Grid: {len(departure_times)} x {len(tof_times)} = {len(departure_times) * len(tof_times)} calculations")
    
    # Pre-calculate Moon positions
    max_time = dep_duration + tof_max + 7200  # Add 2h buffer
    moon_time_step = 1800  # 30 minute intervals for Moon
    moon_times = np.arange(0, max_time, moon_time_step)
    
    print("Pre-calculating Moon positions...")
    moon_positions = []
    for t in moon_times:
        _, _, x, y = moon.get_position_at_time(t)
        moon_positions.append([x, y, 0])
    moon_positions = np.array(moon_positions)
    
    # Initialize result matrices
    dv_total = np.full((len(tof_times), len(departure_times)), np.nan)
    dv_departure = np.full((len(tof_times), len(departure_times)), np.nan)
    dv_arrival = np.full((len(tof_times), len(departure_times)), np.nan)
    
    print("Calculating transfer solutions...")
    total_calcs = len(departure_times) * len(tof_times)
    
    for i, dep_time in enumerate(departure_times):
        # Progress update per row
        progress = (i + 1) / len(departure_times) * 100
        print(f"Row {i+1}/{len(departure_times)} ({progress:.1f}%)")
        
        for j, tof in enumerate(tof_times):
            try:
                # Get LEO position at departure
                _, _, x1, y1 = leo_satellite.get_position_at_time(dep_time)
                R1 = np.array([x1, y1, 0])
                
                # Get Moon position at arrival (interpolate)
                arrival_time = dep_time + tof
                idx = np.searchsorted(moon_times, arrival_time)
                
                if idx == 0:
                    R2 = moon_positions[0]
                elif idx >= len(moon_positions):
                    R2 = moon_positions[-1]
                else:
                    # Linear interpolation
                    t1, t2 = moon_times[idx-1], moon_times[idx]
                    p1, p2 = moon_positions[idx-1], moon_positions[idx]
                    alpha = (arrival_time - t1) / (t2 - t1)
                    R2 = p1 + alpha * (p2 - p1)
                
                # Solve Lambert's problem
                V1, V2, _, _ = lambert.lambert(R1, R2, tof, 'pro', MU_EARTH_KM3_S2)
                
                # Calculate delta-V components
                r1 = np.linalg.norm(R1)
                r2 = np.linalg.norm(R2)
                
                v1_circ = math.sqrt(MU_EARTH_KM3_S2 / r1)
                v2_circ = math.sqrt(MU_EARTH_KM3_S2 / r2)
                
                v1_trans = np.linalg.norm(V1)
                v2_trans = np.linalg.norm(V2)
                
                dv_dep = abs(v1_trans - v1_circ)
                dv_arr = abs(v2_trans - v2_circ)
                dv_tot = dv_dep + dv_arr
                
                # Store if reasonable
                if dv_tot < 15:  # Upper limit
                    dv_total[j, i] = dv_tot
                    dv_departure[j, i] = dv_dep
                    dv_arrival[j, i] = dv_arr
                
            except Exception:
                # Skip failed calculations
                continue
    
    print("Calculations complete! Creating plots...")
    
    # Convert to plotting coordinates
    dep_grid, tof_grid = np.meshgrid(departure_times / 60, tof_times / 3600)
    
    # Create plots
    def create_plot(data_matrix, title, filename):
        """Create and save a porkchop plot"""
        fig, ax = plt.subplots(figsize=(14, 10))
        
        valid_mask = ~np.isnan(data_matrix)
        
        if np.any(valid_mask):
            # Color scale
            vmin = np.nanpercentile(data_matrix, 5)
            vmax = np.nanpercentile(data_matrix, 95)
            
            # Main contour plot
            levels = np.linspace(vmin, vmax, 25)
            contour = ax.contourf(dep_grid, tof_grid, data_matrix, 
                                levels=levels, cmap='plasma', extend='both')
            
            # Contour lines with labels
            contour_lines = ax.contour(dep_grid, tof_grid, data_matrix, 
                                     levels=15, colors='white', alpha=0.7, linewidths=0.8)
            ax.clabel(contour_lines, inline=True, fontsize=10, fmt='%.2f')
            
            # Colorbar
            cbar = plt.colorbar(contour, ax=ax, shrink=0.8)
            cbar.set_label('ΔV (km/s)', rotation=270, labelpad=20, fontsize=12)
            
            # Find optimal solution
            min_idx = np.unravel_index(np.nanargmin(data_matrix), data_matrix.shape)
            opt_dep = dep_grid[min_idx]
            opt_tof = tof_grid[min_idx]
            opt_dv = data_matrix[min_idx]
            
            # Mark optimal point
            ax.plot(opt_dep, opt_tof, 'r*', markersize=20, 
                   markeredgecolor='white', markeredgewidth=2,
                   label=f'Optimal: {opt_dv:.3f} km/s')
            
            # Add LEO orbit reference lines
            leo_period_min = leo_satellite.period / 60
            for n in range(1, 4):
                ax.axvline(n * leo_period_min, color='cyan', linestyle='--', 
                          alpha=0.8, linewidth=1)
                ax.text(n * leo_period_min + 2, ax.get_ylim()[1] * 0.95, 
                       f'{n} orbit', rotation=90, color='cyan', fontsize=10)
            
            # Statistics text box
            stats_text = f'Optimal Solution:\\nDeparture: {opt_dep:.1f} min\\nTOF: {opt_tof:.1f} h\\nΔV: {opt_dv:.3f} km/s'
            ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
                   verticalalignment='top', fontsize=11,
                   bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9))
            
            print(f"Plot {filename}: Optimal ΔV = {opt_dv:.3f} km/s at {opt_dep:.1f} min, {opt_tof:.1f} h")
            
        else:
            ax.text(0.5, 0.5, 'No valid solutions found', 
                   transform=ax.transAxes, ha='center', va='center', fontsize=16)
        
        # Formatting
        ax.set_xlabel('Departure Time (minutes from initial LEO position)', fontsize=12)
        ax.set_ylabel('Time of Flight (hours)', fontsize=12)
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        if 'opt_dv' in locals():
            ax.legend(fontsize=11)
        
        plt.tight_layout()
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved: {filename}")
    
    # Create all three plots
    create_plot(dv_total, 'LEO-to-Moon Transfer: Total ΔV', 'porkchop_total_deltav.png')
    create_plot(dv_departure, 'LEO-to-Moon Transfer: Departure ΔV', 'porkchop_departure_deltav.png')
    create_plot(dv_arrival, 'LEO-to-Moon Transfer: Arrival ΔV', 'porkchop_arrival_deltav.png')
    
    # Summary statistics
    if not np.all(np.isnan(dv_total)):
        min_idx = np.unravel_index(np.nanargmin(dv_total), dv_total.shape)
        
        print(f"\\nFINAL OPTIMAL SOLUTION:")
        print(f"Departure time: {dep_grid[min_idx]:.1f} minutes")
        print(f"Time of flight: {tof_grid[min_idx]:.1f} hours ({tof_grid[min_idx]/24:.1f} days)")
        print(f"Total ΔV: {dv_total[min_idx]:.3f} km/s")
        print(f"  - Departure ΔV: {dv_departure[min_idx]:.3f} km/s")
        print(f"  - Arrival ΔV: {dv_arrival[min_idx]:.3f} km/s")
        
        # Additional analysis
        leo_orbit_fraction = (dep_grid[min_idx] % (leo_satellite.period/60)) / (leo_satellite.period/60)
        print(f"Departure at {leo_orbit_fraction*100:.1f}% through LEO orbit")
    
    return dv_total, dv_departure, dv_arrival, dep_grid, tof_grid

def main():
    """Main function"""
    results = generate_porkchop_plot_with_save()
    print("\\nPorkchop analysis complete! Check the generated PNG files.")
    return results

if __name__ == "__main__":
    results = main()