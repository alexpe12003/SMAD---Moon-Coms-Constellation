"""
Medium Resolution Porkchop Plot
==============================
Optimized version with reasonable resolution for practical analysis.
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import sys
import os

# Add the lambert solver to the path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'Duda_Lambert', 'lambert')))

import lambert
from porkchop import create_leo_satellite, create_moon_orbit

# Physical Constants
MU_EARTH_KM3_S2 = 398600.4418

def medium_resolution_porkchop():
    """Generate a medium resolution porkchop plot"""
    
    print("MEDIUM RESOLUTION PORKCHOP PLOT")
    print("=" * 40)
    
    # Create orbit objects
    leo_satellite = create_leo_satellite()
    moon = create_moon_orbit()
    
    print(f"LEO period: {leo_satellite.period/60:.1f} minutes")
    print(f"Moon period: {moon.period/3600/24:.1f} days")
    
    # Define parameters - optimized for reasonable computation time
    dep_duration = 2 * leo_satellite.period  # 2 LEO orbits
    tof_min = 24 * 3600                      # 1 day minimum
    tof_max = 96 * 3600                      # 4 days maximum
    
    # Grids with moderate resolution
    departure_times = np.linspace(0, dep_duration, 20)  # 20 departure points
    tof_times = np.linspace(tof_min, tof_max, 15)       # 15 time of flight points
    
    print(f"Grid: {len(departure_times)} x {len(tof_times)} = {len(departure_times) * len(tof_times)} calculations")
    
    # Pre-calculate Moon positions for efficiency
    max_time = dep_duration + tof_max + 3600  # Add buffer
    moon_times = np.arange(0, max_time, 3600)  # Hourly Moon positions
    moon_positions = []
    
    print("Pre-calculating Moon positions...")
    for t in moon_times:
        _, _, x, y = moon.get_position_at_time(t)
        moon_positions.append([x, y, 0])
    moon_positions = np.array(moon_positions)
    
    # Initialize result matrix
    dv_matrix = np.full((len(tof_times), len(departure_times)), np.nan)
    
    print("Calculating transfer delta-V...")
    total_calcs = len(departure_times) * len(tof_times)
    
    for i, dep_time in enumerate(departure_times):
        for j, tof in enumerate(tof_times):
            calc_num = i * len(tof_times) + j + 1
            
            if calc_num % 50 == 0:
                print(f"Progress: {calc_num}/{total_calcs} ({calc_num/total_calcs*100:.1f}%)")
            
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
                
                # Calculate delta-V
                r1 = np.linalg.norm(R1)
                r2 = np.linalg.norm(R2)
                
                v1_circ = math.sqrt(MU_EARTH_KM3_S2 / r1)
                v2_circ = math.sqrt(MU_EARTH_KM3_S2 / r2)
                
                v1_trans = np.linalg.norm(V1)
                v2_trans = np.linalg.norm(V2)
                
                dv1 = abs(v1_trans - v1_circ)
                dv2 = abs(v2_trans - v2_circ)
                dv_total = dv1 + dv2
                
                # Store if reasonable
                if dv_total < 12:  # Reasonable upper limit
                    dv_matrix[j, i] = dv_total
                
            except Exception as e:
                # Skip failed calculations
                pass
    
    print("Calculation complete! Creating plot...")
    
    # Create plot
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Convert to plotting units
    dep_grid, tof_grid = np.meshgrid(departure_times / 60, tof_times / 3600)
    
    # Plot only if we have valid data
    valid_mask = ~np.isnan(dv_matrix)
    
    if np.any(valid_mask):
        # Create contour plot
        min_dv = np.nanmin(dv_matrix)
        max_dv = np.nanmax(dv_matrix)
        
        levels = np.linspace(min_dv, max_dv, 20)
        contour = ax.contourf(dep_grid, tof_grid, dv_matrix, 
                            levels=levels, cmap='plasma')
        
        # Add contour lines with labels
        contour_lines = ax.contour(dep_grid, tof_grid, dv_matrix, 
                                 levels=10, colors='white', alpha=0.6)
        ax.clabel(contour_lines, inline=True, fontsize=8, fmt='%.2f')
        
        # Colorbar
        cbar = plt.colorbar(contour, ax=ax)
        cbar.set_label('Total ΔV (km/s)', rotation=270, labelpad=20)
        
        # Find and mark minimum
        min_idx = np.unravel_index(np.nanargmin(dv_matrix), dv_matrix.shape)
        min_dep = dep_grid[min_idx]
        min_tof = tof_grid[min_idx]
        min_dv_val = dv_matrix[min_idx]
        
        ax.plot(min_dep, min_tof, 'r*', markersize=20, 
               label=f'Optimal: {min_dv_val:.3f} km/s')
        
        # Add LEO orbit reference lines
        leo_period_min = leo_satellite.period / 60
        ax.axvline(leo_period_min, color='red', linestyle='--', alpha=0.7, 
                  label=f'1 LEO orbit ({leo_period_min:.1f} min)')
        ax.axvline(2 * leo_period_min, color='red', linestyle='--', alpha=0.7)
        
        # Statistics text
        ax.text(0.02, 0.98, f'Optimal Solution:\nDeparture: {min_dep:.1f} min\nTOF: {min_tof:.1f} h\nΔV: {min_dv_val:.3f} km/s', 
               transform=ax.transAxes, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
        
        print(f"\nOPTIMAL SOLUTION FOUND:")
        print(f"Departure time: {min_dep:.1f} minutes")
        print(f"Time of flight: {min_tof:.1f} hours ({min_tof/24:.1f} days)")
        print(f"Total ΔV: {min_dv_val:.3f} km/s")
        
    else:
        ax.text(0.5, 0.5, 'No valid transfer solutions found', 
               transform=ax.transAxes, ha='center', va='center')
        print("No valid solutions found!")
    
    # Formatting
    ax.set_xlabel('Departure Time (minutes from initial LEO position)')
    ax.set_ylabel('Time of Flight (hours)')
    ax.set_title('LEO-to-Moon Transfer Porkchop Plot\n(Medium Resolution)')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    plt.tight_layout()
    plt.show()
    
    return dv_matrix, dep_grid, tof_grid

def main():
    """Main function"""
    result = medium_resolution_porkchop()
    print("\nPorkchop plot generation complete!")
    return result

if __name__ == "__main__":
    result = main()