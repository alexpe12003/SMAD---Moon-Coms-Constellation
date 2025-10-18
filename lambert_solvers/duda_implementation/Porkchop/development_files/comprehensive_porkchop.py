"""
Comprehensive Porkchop Plot Generator
====================================
Enhanced version for generating high-quality porkchop plots for LEO-to-Moon transfers.
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

class PorkchopGenerator:
    """Enhanced porkchop plot generator"""
    
    def __init__(self):
        self.leo_satellite = create_leo_satellite()
        self.moon = create_moon_orbit()
        
        print("PORKCHOP GENERATOR INITIALIZED")
        print("=" * 40)
        print(f"LEO period: {self.leo_satellite.period/60:.1f} minutes")
        print(f"Moon period: {self.moon.period/3600/24:.1f} days")
    
    def extend_moon_orbit(self, duration_days=7):
        """Generate extended Moon orbit data"""
        total_time = duration_days * 24 * 3600  # Convert to seconds
        time_step = 3600  # 1 hour steps for Moon (it moves slowly)
        
        time_points = np.arange(0, total_time + time_step, time_step)
        
        positions = []
        for t in time_points:
            _, _, x, y = self.moon.get_position_at_time(t)
            positions.append([x, y, 0])  # Z=0 for 2D
        
        return np.array(time_points), np.array(positions)
    
    def get_moon_position_at_time(self, time_s, moon_times, moon_positions):
        """Interpolate Moon position at any time"""
        if time_s <= moon_times[0]:
            return moon_positions[0]
        elif time_s >= moon_times[-1]:
            return moon_positions[-1]
        else:
            # Linear interpolation
            idx = np.searchsorted(moon_times, time_s)
            if idx == 0:
                return moon_positions[0]
            
            t1, t2 = moon_times[idx-1], moon_times[idx]
            p1, p2 = moon_positions[idx-1], moon_positions[idx]
            
            # Interpolation factor
            alpha = (time_s - t1) / (t2 - t1)
            return p1 + alpha * (p2 - p1)
    
    def calculate_transfer_dv(self, R1, R2, tof):
        """Calculate delta-V for transfer"""
        try:
            V1, V2, _, _ = lambert.lambert(R1, R2, tof, 'pro', MU_EARTH_KM3_S2)
            
            # Circular velocities
            r1 = np.linalg.norm(R1)
            r2 = np.linalg.norm(R2)
            
            v1_circ = math.sqrt(MU_EARTH_KM3_S2 / r1)
            v2_circ = math.sqrt(MU_EARTH_KM3_S2 / r2)
            
            # Transfer velocities
            v1_trans = np.linalg.norm(V1)
            v2_trans = np.linalg.norm(V2)
            
            # Delta-V components
            dv1 = abs(v1_trans - v1_circ)
            dv2 = abs(v2_trans - v2_circ)
            
            return dv1 + dv2, dv1, dv2, V1, V2
            
        except Exception as e:
            return np.nan, np.nan, np.nan, None, None
    
    def generate_porkchop_data(self, dep_duration_orbits=4, tof_min_hours=12, 
                             tof_max_hours=120, dep_resolution_min=5, tof_resolution_hours=2):
        """Generate comprehensive porkchop data"""
        
        print(f"\nGenerating porkchop data...")
        print(f"Departure: {dep_duration_orbits} LEO orbits")
        print(f"Time of flight: {tof_min_hours} to {tof_max_hours} hours")
        
        # Generate extended Moon data
        moon_sim_days = max(7, tof_max_hours/24 + 2)  # Ensure we have enough Moon data
        moon_times, moon_positions = self.extend_moon_orbit(moon_sim_days)
        print(f"Moon simulation: {moon_sim_days} days")
        
        # Define grids
        dep_max_time = dep_duration_orbits * self.leo_satellite.period
        departure_times = np.arange(0, dep_max_time, dep_resolution_min * 60)
        tof_times = np.arange(tof_min_hours * 3600, (tof_max_hours + tof_resolution_hours) * 3600, 
                             tof_resolution_hours * 3600)
        
        print(f"Grid size: {len(departure_times)} x {len(tof_times)} = {len(departure_times) * len(tof_times)} calculations")
        
        # Initialize matrices
        dv_total = np.full((len(tof_times), len(departure_times)), np.nan)
        dv_departure = np.full((len(tof_times), len(departure_times)), np.nan)
        dv_arrival = np.full((len(tof_times), len(departure_times)), np.nan)
        
        # Calculate delta-V for each combination
        total_calcs = len(departure_times) * len(tof_times)
        calc_count = 0
        
        for i, dep_time in enumerate(departure_times):
            for j, tof in enumerate(tof_times):
                calc_count += 1
                
                # Progress update
                if calc_count % 200 == 0:
                    print(f"Progress: {calc_count/total_calcs*100:.1f}%")
                
                # Get LEO position at departure
                _, _, x1, y1 = self.leo_satellite.get_position_at_time(dep_time)
                R1 = np.array([x1, y1, 0])
                
                # Get Moon position at arrival
                arrival_time = dep_time + tof
                R2 = self.get_moon_position_at_time(arrival_time, moon_times, moon_positions)
                
                # Calculate delta-V
                dv_tot, dv_dep, dv_arr, V1, V2 = self.calculate_transfer_dv(R1, R2, tof)
                
                # Apply reasonable constraints
                if not np.isnan(dv_tot) and dv_tot < 15:  # Reasonable upper limit
                    dv_total[j, i] = dv_tot
                    dv_departure[j, i] = dv_dep
                    dv_arrival[j, i] = dv_arr
        
        print("Data generation complete!")
        
        # Create coordinate grids for plotting
        dep_grid, tof_grid = np.meshgrid(departure_times / 60, tof_times / 3600)  # Convert to minutes and hours
        
        return dv_total, dv_departure, dv_arrival, dep_grid, tof_grid
    
    def plot_porkchop(self, dv_matrix, dep_grid, tof_grid, title, vmin=None, vmax=None):
        """Create a porkchop plot"""
        
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Filter valid data
        valid_mask = ~np.isnan(dv_matrix)
        
        if np.any(valid_mask):
            # Set color scale limits
            if vmin is None:
                vmin = np.nanpercentile(dv_matrix, 5)
            if vmax is None:
                vmax = np.nanpercentile(dv_matrix, 95)
            
            # Create contour plot
            levels = np.linspace(vmin, vmax, 25)
            contour = ax.contourf(dep_grid, tof_grid, dv_matrix, 
                                levels=levels, cmap='viridis', extend='both')
            
            # Add contour lines
            contour_lines = ax.contour(dep_grid, tof_grid, dv_matrix, 
                                     levels=levels, colors='white', alpha=0.4, linewidths=0.5)
            
            # Add specific contour labels for key values
            key_levels = [3, 4, 5, 6, 7, 8]
            key_contours = ax.contour(dep_grid, tof_grid, dv_matrix, 
                                    levels=key_levels, colors='black', linewidths=1)
            ax.clabel(key_contours, inline=True, fontsize=8, fmt='%.1f')
            
            # Colorbar
            cbar = plt.colorbar(contour, ax=ax)
            cbar.set_label('ΔV (km/s)', rotation=270, labelpad=20)
            
            # Find and mark minimum
            min_idx = np.unravel_index(np.nanargmin(dv_matrix), dv_matrix.shape)
            min_dep = dep_grid[min_idx]
            min_tof = tof_grid[min_idx]
            min_dv = dv_matrix[min_idx]
            
            ax.plot(min_dep, min_tof, 'r*', markersize=15, 
                   label=f'Minimum: {min_dv:.3f} km/s')
            
            # Add LEO orbit markers
            leo_period_min = self.leo_satellite.period / 60
            for i in range(1, int(dep_grid.max() / leo_period_min) + 1):
                ax.axvline(i * leo_period_min, color='red', linestyle='--', alpha=0.5)
            
            ax.text(0.02, 0.98, f'Min ΔV: {min_dv:.3f} km/s\nAt t={min_dep:.1f} min\nTOF={min_tof:.1f} h', 
                   transform=ax.transAxes, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
            
        else:
            ax.text(0.5, 0.5, 'No valid solutions found', 
                   transform=ax.transAxes, ha='center', va='center')
        
        ax.set_xlabel('Departure Time (minutes from initial LEO position)')
        ax.set_ylabel('Time of Flight (hours)')
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        plt.tight_layout()
        return fig, ax
    
    def generate_full_analysis(self):
        """Generate complete porkchop analysis"""
        
        print("GENERATING COMPREHENSIVE PORKCHOP ANALYSIS")
        print("=" * 50)
        
        # Generate data with good resolution
        dv_total, dv_dep, dv_arr, dep_grid, tof_grid = self.generate_porkchop_data(
            dep_duration_orbits=4,     # 4 LEO orbits
            tof_min_hours=12,          # Minimum 12 hours
            tof_max_hours=120,         # Maximum 5 days
            dep_resolution_min=3,      # Every 3 minutes
            tof_resolution_hours=1     # Every hour
        )
        
        # Create plots
        print("\nCreating plots...")
        
        # Total Delta-V plot
        fig1, ax1 = self.plot_porkchop(dv_total, dep_grid, tof_grid, 
                                      'LEO-to-Moon Transfer: Total ΔV')
        
        # Departure Delta-V plot
        fig2, ax2 = self.plot_porkchop(dv_dep, dep_grid, tof_grid, 
                                      'LEO-to-Moon Transfer: Departure ΔV')
        
        # Arrival Delta-V plot
        fig3, ax3 = self.plot_porkchop(dv_arr, dep_grid, tof_grid, 
                                      'LEO-to-Moon Transfer: Arrival ΔV')
        
        plt.show()
        
        # Find and report optimal solution
        if not np.all(np.isnan(dv_total)):
            min_idx = np.unravel_index(np.nanargmin(dv_total), dv_total.shape)
            optimal_dep = dep_grid[min_idx]
            optimal_tof = tof_grid[min_idx]
            optimal_dv_total = dv_total[min_idx]
            optimal_dv_dep = dv_dep[min_idx]
            optimal_dv_arr = dv_arr[min_idx]
            
            print(f"\nOPTIMAL TRANSFER SOLUTION:")
            print(f"Departure time: {optimal_dep:.1f} minutes ({optimal_dep/self.leo_satellite.period*60:.1f}% of LEO orbit)")
            print(f"Time of flight: {optimal_tof:.1f} hours ({optimal_tof/24:.1f} days)")
            print(f"Total ΔV: {optimal_dv_total:.3f} km/s")
            print(f"  - Departure ΔV: {optimal_dv_dep:.3f} km/s")
            print(f"  - Arrival ΔV: {optimal_dv_arr:.3f} km/s")
        
        return dv_total, dv_dep, dv_arr, dep_grid, tof_grid

def main():
    """Main function"""
    generator = PorkchopGenerator()
    results = generator.generate_full_analysis()
    
    print("\nPorkchop analysis complete!")
    return results

if __name__ == "__main__":
    results = main()