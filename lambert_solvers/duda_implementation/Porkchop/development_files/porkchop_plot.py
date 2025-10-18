"""
Porkchop Plot Generator for LEO-to-Moon Transfer
===============================================
This script generates a porkchop plot showing the delta-V requirements
for transfers from LEO orbit to Moon orbit for various departure times
and time of flight durations.

The porkchop plot displays:
- X-axis: Departure time (position in LEO orbit)
- Y-axis: Time of flight 
- Contours: Total delta-V required for the transfer

Author: Generated for SMAD Moon Communications Constellation
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import sys
import os

# Add the lambert solver to the path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'Duda_Lambert', 'lambert')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'Duda_Lambert', 'Auxiliary')))

import lambert  # Import the lambert module
from porkchop import create_leo_satellite, create_moon_orbit

# Physical Constants
MU_EARTH_KM3_S2 = 398600.4418  # Earth's gravitational parameter (km³/s²)

class PorkchopAnalyzer:
    """Class to generate porkchop plots for LEO-to-Moon transfers"""
    
    def __init__(self):
        """Initialize the porkchop analyzer"""
        self.leo_satellite = create_leo_satellite()
        self.moon = create_moon_orbit()
        
        print("PORKCHOP PLOT ANALYZER INITIALIZED")
        print("=" * 50)
        print(f"LEO orbital period: {self.leo_satellite.period/60:.1f} minutes")
        print(f"Moon orbital period: {self.moon.period/3600/24:.1f} days")
        
    def extend_moon_orbit(self, duration_hours=168):  # 1 week default
        """
        Generate extended Moon orbit data to cover longer time periods
        
        Parameters:
        - duration_hours: Duration to simulate in hours
        
        Returns:
        - Dictionary with extended moon orbit data
        """
        print(f"\nGenerating extended Moon orbit data for {duration_hours} hours...")
        
        time_interval_s = 60  # 1 minute intervals
        total_time_s = duration_hours * 3600
        time_points = np.arange(0, total_time_s + time_interval_s, time_interval_s)
        
        moon_data = {
            'time_s': time_points,
            'time_min': time_points / 60,
            'time_h': time_points / 3600,
            'theta_rad': [],
            'theta_deg': [],
            'x_km': [],
            'y_km': []
        }
        
        for t in time_points:
            theta_rad, theta_deg, x_km, y_km = self.moon.get_position_at_time(t)
            moon_data['theta_rad'].append(theta_rad)
            moon_data['theta_deg'].append(theta_deg)
            moon_data['x_km'].append(x_km)
            moon_data['y_km'].append(y_km)
        
        # Convert to numpy arrays
        for key in ['theta_rad', 'theta_deg', 'x_km', 'y_km']:
            moon_data[key] = np.array(moon_data[key])
        
        print(f"Generated {len(time_points)} Moon position data points")
        return moon_data
    
    def get_leo_position(self, departure_time_s):
        """
        Get LEO satellite position at departure time
        
        Parameters:
        - departure_time_s: Departure time in seconds
        
        Returns:
        - R1: Position vector [x, y, z] in km
        """
        theta_rad, theta_deg, x_km, y_km = self.leo_satellite.get_position_at_time(departure_time_s)
        return np.array([x_km, y_km, 0])  # Z=0 for 2D coplanar orbits
    
    def get_moon_position(self, time_s, moon_data):
        """
        Get Moon position at specified time using interpolation
        
        Parameters:
        - time_s: Time in seconds
        - moon_data: Extended moon orbit data
        
        Returns:
        - R2: Position vector [x, y, z] in km
        """
        # Find closest time index or interpolate
        time_idx = np.searchsorted(moon_data['time_s'], time_s)
        
        if time_idx >= len(moon_data['time_s']):
            time_idx = len(moon_data['time_s']) - 1
        
        x_km = moon_data['x_km'][time_idx]
        y_km = moon_data['y_km'][time_idx]
        
        return np.array([x_km, y_km, 0])  # Z=0 for 2D coplanar orbits
    
    def calculate_transfer_delta_v(self, R1, R2, tof_s):
        """
        Calculate delta-V for transfer using Lambert's problem
        
        Parameters:
        - R1: Initial position vector (km)
        - R2: Final position vector (km)
        - tof_s: Time of flight in seconds
        
        Returns:
        - delta_v_total: Total delta-V (km/s)
        - v1_transfer: Transfer velocity at departure (km/s)
        - v2_transfer: Transfer velocity at arrival (km/s)
        """
        try:
            # Solve Lambert's problem
            V1_transfer, V2_transfer, theta1, theta2 = lambert.lambert(R1, R2, tof_s, 'pro', MU_EARTH_KM3_S2)
            
            # Calculate circular orbital velocities
            r1 = np.linalg.norm(R1)
            r2 = np.linalg.norm(R2)
            
            v1_circular = math.sqrt(MU_EARTH_KM3_S2 / r1)  # LEO circular velocity
            v2_circular = math.sqrt(MU_EARTH_KM3_S2 / r2)  # Moon circular velocity
            
            # Calculate delta-V requirements
            # Delta-V at departure (LEO)
            v1_transfer_mag = np.linalg.norm(V1_transfer)
            delta_v1 = abs(v1_transfer_mag - v1_circular)
            
            # Delta-V at arrival (Moon)
            v2_transfer_mag = np.linalg.norm(V2_transfer)
            delta_v2 = abs(v2_transfer_mag - v2_circular)
            
            # Total delta-V
            delta_v_total = delta_v1 + delta_v2
            
            return delta_v_total, V1_transfer, V2_transfer, delta_v1, delta_v2
            
        except Exception as e:
            print(f"Error in Lambert solver: {e}")
            return np.nan, None, None, np.nan, np.nan
    
    def generate_porkchop_data(self, departure_times_min, tof_range_hours, moon_data):
        """
        Generate porkchop plot data matrix
        
        Parameters:
        - departure_times_min: Array of departure times in minutes
        - tof_range_hours: Array of time of flight values in hours
        - moon_data: Extended moon orbit data
        
        Returns:
        - delta_v_matrix: 2D array of delta-V values
        - departure_grid: 2D grid of departure times
        - tof_grid: 2D grid of time of flight values
        """
        print(f"\nGenerating porkchop data...")
        print(f"Departure times: {len(departure_times_min)} points")
        print(f"Time of flight range: {len(tof_range_hours)} points")
        
        # Initialize matrices
        delta_v_matrix = np.zeros((len(tof_range_hours), len(departure_times_min)))
        delta_v1_matrix = np.zeros((len(tof_range_hours), len(departure_times_min)))
        delta_v2_matrix = np.zeros((len(tof_range_hours), len(departure_times_min)))
        
        total_calculations = len(departure_times_min) * len(tof_range_hours)
        calculation_count = 0
        
        for i, departure_time_min in enumerate(departure_times_min):
            for j, tof_hours in enumerate(tof_range_hours):
                calculation_count += 1
                
                # Progress indicator
                if calculation_count % 100 == 0:
                    progress = (calculation_count / total_calculations) * 100
                    print(f"Progress: {progress:.1f}% ({calculation_count}/{total_calculations})")
                
                # Convert to seconds
                departure_time_s = departure_time_min * 60
                tof_s = tof_hours * 3600
                arrival_time_s = departure_time_s + tof_s
                
                # Get positions
                R1 = self.get_leo_position(departure_time_s)
                R2 = self.get_moon_position(arrival_time_s, moon_data)
                
                # Calculate delta-V
                delta_v_total, V1, V2, delta_v1, delta_v2 = self.calculate_transfer_delta_v(R1, R2, tof_s)
                
                # Store results
                delta_v_matrix[j, i] = delta_v_total
                delta_v1_matrix[j, i] = delta_v1
                delta_v2_matrix[j, i] = delta_v2
        
        print("Porkchop data generation complete!")
        
        # Create coordinate grids for plotting
        departure_grid, tof_grid = np.meshgrid(departure_times_min, tof_range_hours)
        
        return delta_v_matrix, delta_v1_matrix, delta_v2_matrix, departure_grid, tof_grid
    
    def plot_porkchop(self, delta_v_matrix, departure_grid, tof_grid, 
                     title="LEO-to-Moon Transfer Porkchop Plot"):
        """
        Create the porkchop plot
        
        Parameters:
        - delta_v_matrix: 2D array of delta-V values
        - departure_grid: 2D grid of departure times
        - tof_grid: 2D grid of time of flight values
        - title: Plot title
        """
        plt.figure(figsize=(12, 8))
        
        # Create contour plot
        # Filter out NaN values for better visualization
        valid_mask = ~np.isnan(delta_v_matrix)
        if np.any(valid_mask):
            min_dv = np.nanmin(delta_v_matrix)
            max_dv = np.nanmax(delta_v_matrix)
            
            # Create contour levels
            contour_levels = np.linspace(min_dv, max_dv, 20)
            
            # Main contour plot
            contour = plt.contourf(departure_grid, tof_grid, delta_v_matrix, 
                                 levels=contour_levels, cmap='viridis')
            
            # Add contour lines
            contour_lines = plt.contour(departure_grid, tof_grid, delta_v_matrix, 
                                      levels=contour_levels, colors='white', alpha=0.5, linewidths=0.5)
            
            # Add colorbar
            cbar = plt.colorbar(contour)
            cbar.set_label('Total ΔV (km/s)', rotation=270, labelpad=20)
            
            # Labels and formatting
            plt.xlabel('Departure Time (minutes from initial LEO position)')
            plt.ylabel('Time of Flight (hours)')
            plt.title(title)
            plt.grid(True, alpha=0.3)
            
            # Add some reference lines
            plt.axvline(x=self.leo_satellite.period/60, color='red', linestyle='--', alpha=0.7, 
                       label=f'One LEO orbit ({self.leo_satellite.period/60:.1f} min)')
            
            plt.legend()
            
            # Statistics
            plt.text(0.02, 0.98, f'Min ΔV: {min_dv:.3f} km/s\nMax ΔV: {max_dv:.3f} km/s', 
                    transform=plt.gca().transAxes, verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        else:
            plt.text(0.5, 0.5, 'No valid transfer solutions found', 
                    transform=plt.gca().transAxes, ha='center', va='center')
        
        plt.tight_layout()
        plt.show()
    
    def generate_full_porkchop_plot(self, 
                                   departure_duration_hours=6,  # Cover multiple LEO orbits
                                   tof_min_hours=10,            # Minimum time of flight
                                   tof_max_hours=120,           # Maximum time of flight (5 days)
                                   moon_sim_duration_hours=168, # 1 week of Moon data
                                   resolution_factor=1):         # 1=full resolution, 2=half, etc.
        """
        Generate a complete porkchop plot with specified parameters
        
        Parameters:
        - departure_duration_hours: How many hours of departure times to cover
        - tof_min_hours: Minimum time of flight
        - tof_max_hours: Maximum time of flight  
        - moon_sim_duration_hours: Duration to simulate Moon orbit
        - resolution_factor: Factor to reduce resolution (for faster calculation)
        """
        print(f"\nGENERATING PORKCHOP PLOT")
        print("=" * 50)
        print(f"Departure duration: {departure_duration_hours} hours")
        print(f"Time of flight range: {tof_min_hours} to {tof_max_hours} hours")
        print(f"Resolution factor: {resolution_factor}")
        
        # Generate extended Moon orbit data
        moon_data = self.extend_moon_orbit(moon_sim_duration_hours)
        
        # Define departure times and time of flight ranges
        departure_times_min = np.arange(0, departure_duration_hours * 60, 
                                      10 * resolution_factor)  # Every 10 minutes
        tof_range_hours = np.arange(tof_min_hours, tof_max_hours + 1, 
                                  2 * resolution_factor)    # Every 2 hours
        
        print(f"Departure times: {len(departure_times_min)} points")
        print(f"Time of flight: {len(tof_range_hours)} points")
        print(f"Total calculations: {len(departure_times_min) * len(tof_range_hours)}")
        
        # Generate porkchop data
        delta_v_matrix, delta_v1_matrix, delta_v2_matrix, departure_grid, tof_grid = \
            self.generate_porkchop_data(departure_times_min, tof_range_hours, moon_data)
        
        # Create plots
        # Main porkchop plot
        self.plot_porkchop(delta_v_matrix, departure_grid, tof_grid,
                          "LEO-to-Moon Transfer Porkchop Plot (Total ΔV)")
        
        # Departure delta-V plot
        self.plot_porkchop(delta_v1_matrix, departure_grid, tof_grid,
                          "LEO-to-Moon Transfer Porkchop Plot (Departure ΔV)")
        
        # Arrival delta-V plot
        self.plot_porkchop(delta_v2_matrix, departure_grid, tof_grid,
                          "LEO-to-Moon Transfer Porkchop Plot (Arrival ΔV)")
        
        return delta_v_matrix, delta_v1_matrix, delta_v2_matrix, departure_grid, tof_grid

def main():
    """Main function to generate porkchop plots"""
    
    print("LEO-TO-MOON TRANSFER PORKCHOP PLOT GENERATOR")
    print("=" * 60)
    
    # Create analyzer
    analyzer = PorkchopAnalyzer()
    
    # Generate porkchop plot with moderate resolution for testing
    print(f"\nGenerating porkchop plot with moderate resolution...")
    
    results = analyzer.generate_full_porkchop_plot(
        departure_duration_hours=6,    # Cover ~4 LEO orbits
        tof_min_hours=10,             # Minimum 10 hours
        tof_max_hours=120,            # Maximum 5 days
        moon_sim_duration_hours=168,  # 1 week of Moon simulation
        resolution_factor=2           # Half resolution for faster calculation
    )
    
    delta_v_matrix, delta_v1_matrix, delta_v2_matrix, departure_grid, tof_grid = results
    
    # Find minimum delta-V solution
    valid_mask = ~np.isnan(delta_v_matrix)
    if np.any(valid_mask):
        min_idx = np.unravel_index(np.nanargmin(delta_v_matrix), delta_v_matrix.shape)
        min_departure = departure_grid[min_idx]
        min_tof = tof_grid[min_idx]
        min_dv = delta_v_matrix[min_idx]
        
        print(f"\nOPTIMAL TRANSFER SOLUTION:")
        print(f"Departure time: {min_departure:.1f} minutes")
        print(f"Time of flight: {min_tof:.1f} hours")
        print(f"Total ΔV: {min_dv:.3f} km/s")
        print(f"Departure ΔV: {delta_v1_matrix[min_idx]:.3f} km/s")
        print(f"Arrival ΔV: {delta_v2_matrix[min_idx]:.3f} km/s")
    
    print(f"\nPorkchop plot generation complete!")
    
    return results

if __name__ == "__main__":
    results = main()