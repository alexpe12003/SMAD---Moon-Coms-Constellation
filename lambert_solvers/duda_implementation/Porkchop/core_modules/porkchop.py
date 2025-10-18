"""
Circular Orbit Position and Time Calculator
==================================================
This script calculates the position and time for:
1. A satellite in a circular 200km LEO orbit around Earth
2. The Moon in a circular orbit around Earth

Both orbits are assumed to be coplanar and 2D.
Position is defined by the angle θ, with θ=0 at initial position.
Time intervals: 1 minute

Author: Generated for SMAD Moon Communications Constellation
"""

import numpy as np
import matplotlib.pyplot as plt
import math

# Physical Constants
MU_EARTH_KM3_S2 = 398600.4418  # Earth's gravitational parameter (km³/s²)
EARTH_RADIUS_KM = 6378  # Earth's radius (km)
MOON_DISTANCE_KM = 384400  # Average Earth-Moon distance (km)

class CircularOrbit:
    """Class to represent a circular orbit and calculate positions over time"""
    
    def __init__(self, radius_km, mu_km3_s2, name="Orbit"):
        """
        Initialize circular orbit
        
        Parameters:
        - radius_km: Orbital radius in km
        - mu_km3_s2: Gravitational parameter in km³/s²
        - name: Name identifier for the orbit
        """
        self.radius = radius_km
        self.mu = mu_km3_s2
        self.name = name
        
        # Calculate orbital parameters
        self.velocity = math.sqrt(mu_km3_s2 / radius_km)  # Circular orbital velocity
        self.period = 2 * math.pi * math.sqrt(radius_km**3 / mu_km3_s2)  # Orbital period
        self.angular_velocity = 2 * math.pi / self.period  # Angular velocity (rad/s)
        
        print(f"\n{self.name} Orbital Parameters:")
        print(f"  Radius: {self.radius:.0f} km")
        print(f"  Velocity: {self.velocity:.3f} km/s")
        print(f"  Period: {self.period/3600:.2f} hours ({self.period/60:.1f} minutes)")
        print(f"  Angular velocity: {self.angular_velocity:.6f} rad/s")
        print(f"  Angular velocity: {math.degrees(self.angular_velocity):.6f} deg/s")
    
    def get_position_at_time(self, time_s):
        """
        Get position at given time
        
        Parameters:
        - time_s: Time in seconds from initial position
        
        Returns:
        - theta_rad: Angle in radians
        - theta_deg: Angle in degrees
        - x_km: X coordinate in km
        - y_km: Y coordinate in km
        """
        theta_rad = self.angular_velocity * time_s
        theta_deg = math.degrees(theta_rad)
        
        # Position in Cartesian coordinates (assuming θ=0 at (radius, 0))
        x_km = self.radius * math.cos(theta_rad)
        y_km = self.radius * math.sin(theta_rad)
        
        return theta_rad, theta_deg, x_km, y_km
    
    def generate_orbit_data(self, time_interval_s=60, num_orbits=1):
        """
        Generate orbit data over specified time period
        
        Parameters:
        - time_interval_s: Time interval between data points (seconds)
        - num_orbits: Number of complete orbits to calculate
        
        Returns:
        - Dictionary containing time and position arrays
        """
        total_time = num_orbits * self.period
        time_points = np.arange(0, total_time + time_interval_s, time_interval_s)
        
        data = {
            'time_s': time_points,
            'time_min': time_points / 60,
            'time_h': time_points / 3600,
            'theta_rad': [],
            'theta_deg': [],
            'x_km': [],
            'y_km': []
        }
        
        for t in time_points:
            theta_rad, theta_deg, x_km, y_km = self.get_position_at_time(t)
            data['theta_rad'].append(theta_rad)
            data['theta_deg'].append(theta_deg)
            data['x_km'].append(x_km)
            data['y_km'].append(y_km)
        
        # Convert to numpy arrays
        for key in ['theta_rad', 'theta_deg', 'x_km', 'y_km']:
            data[key] = np.array(data[key])
        
        return data

def create_leo_satellite():
    """Create LEO satellite orbit (200km altitude)"""
    leo_radius = EARTH_RADIUS_KM + 200  # 200km altitude
    return CircularOrbit(leo_radius, MU_EARTH_KM3_S2, "LEO Satellite (200km)")

def create_moon_orbit():
    """Create Moon orbit (average Earth-Moon distance)"""
    return CircularOrbit(MOON_DISTANCE_KM, MU_EARTH_KM3_S2, "Moon")

def display_orbit_comparison(leo_data, moon_data, num_points=10):
    """Display comparison of orbit data for first few time points"""
    
    print(f"\n{'='*80}")
    print("ORBIT COMPARISON - First {num_points} Data Points (1-minute intervals)")
    print(f"{'='*80}")
    
    print(f"\n{'Time (min)':<10} {'LEO θ (deg)':<12} {'LEO X (km)':<12} {'LEO Y (km)':<12} {'Moon θ (deg)':<14} {'Moon X (km)':<14} {'Moon Y (km)':<14}")
    print("-" * 80)
    
    for i in range(min(num_points, len(leo_data['time_min']))):
        time_min = leo_data['time_min'][i]
        leo_theta = leo_data['theta_deg'][i] % 360  # Normalize to 0-360
        leo_x = leo_data['x_km'][i]
        leo_y = leo_data['y_km'][i]
        moon_theta = moon_data['theta_deg'][i] % 360  # Normalize to 0-360
        moon_x = moon_data['x_km'][i]
        moon_y = moon_data['y_km'][i]
        
        print(f"{time_min:<10.1f} {leo_theta:<12.1f} {leo_x:<12.0f} {leo_y:<12.0f} {moon_theta:<14.3f} {moon_x:<14.0f} {moon_y:<14.0f}")

def plot_orbits(leo_data, moon_data, time_hours=6):
    """
    Plot the orbits for visualization
    
    Parameters:
    - leo_data: LEO satellite orbit data
    - moon_data: Moon orbit data  
    - time_hours: Time duration to plot (hours)
    """
    # Filter data for specified time duration
    time_mask_leo = leo_data['time_h'] <= time_hours
    time_mask_moon = moon_data['time_h'] <= time_hours
    
    # Calculate actual data ranges
    leo_max_time = np.max(leo_data['time_h'])
    moon_max_time = np.max(moon_data['time_h'])
    leo_plot_time = min(time_hours, leo_max_time)
    moon_plot_time = min(time_hours, moon_max_time)
    
    print(f"📊 Plot Info:")
    print(f"   LEO data available: {leo_max_time:.1f} hours, plotting: {leo_plot_time:.1f} hours")
    print(f"   Moon data available: {moon_max_time:.1f} hours, plotting: {moon_plot_time:.1f} hours")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Plot LEO satellite orbit
    ax1.plot(leo_data['x_km'][time_mask_leo], leo_data['y_km'][time_mask_leo], 'b-', linewidth=2, label='LEO Satellite')
    ax1.plot(leo_data['x_km'][0], leo_data['y_km'][0], 'ro', markersize=8, label='Initial Position')
    ax1.add_patch(plt.Circle((0, 0), EARTH_RADIUS_KM, color='green', alpha=0.3, label='Earth'))
    ax1.set_xlabel('X Position (km)')
    ax1.set_ylabel('Y Position (km)')
    ax1.set_title(f'LEO Satellite Orbit (200km altitude)\n{leo_plot_time:.1f} hours of motion')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.axis('equal')
    
    # Plot Moon orbit (scaled view)
    ax2.plot(moon_data['x_km'][time_mask_moon], moon_data['y_km'][time_mask_moon], 'gray', linewidth=2, label='Moon Orbit')
    ax2.plot(moon_data['x_km'][0], moon_data['y_km'][0], 'yo', markersize=8, label='Initial Position')
    ax2.plot(0, 0, 'go', markersize=10, label='Earth')
    ax2.set_xlabel('X Position (km)')
    ax2.set_ylabel('Y Position (km)')
    ax2.set_title(f'Moon Orbit\n{moon_plot_time:.1f} hours of motion')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    ax2.axis('equal')
    
    plt.tight_layout()
    plt.show()

def plot_angular_positions(leo_data, moon_data, time_hours=24):
    """Plot angular positions vs time"""
    
    # Filter data for specified time duration
    time_mask_leo = leo_data['time_h'] <= time_hours
    time_mask_moon = moon_data['time_h'] <= time_hours
    
    # Calculate actual data ranges
    leo_max_time = np.max(leo_data['time_h'])
    moon_max_time = np.max(moon_data['time_h'])
    leo_plot_time = min(time_hours, leo_max_time)
    moon_plot_time = min(time_hours, moon_max_time)
    
    print(f"📈 Angular Plot Info:")
    print(f"   Plotting time range: {time_hours:.1f} hours")
    print(f"   LEO data: {leo_plot_time:.1f} hours, Moon data: {moon_plot_time:.1f} hours")
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
    
    # LEO satellite angular position
    ax1.plot(leo_data['time_h'][time_mask_leo], leo_data['theta_deg'][time_mask_leo] % 360, 'b-', linewidth=2)
    ax1.set_xlabel('Time (hours)')
    ax1.set_ylabel('Angular Position θ (degrees)')
    ax1.set_title(f'LEO Satellite Angular Position vs Time ({leo_plot_time:.1f} hours)')
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, 360)
    
    # Moon angular position
    ax2.plot(moon_data['time_h'][time_mask_moon], moon_data['theta_deg'][time_mask_moon] % 360, 'gray', linewidth=2)
    ax2.set_xlabel('Time (hours)')
    ax2.set_ylabel('Angular Position θ (degrees)')
    ax2.set_title(f'Moon Angular Position vs Time ({moon_plot_time:.1f} hours)')
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(0, 360)
    
    plt.tight_layout()
    plt.show()

def main():
    """Main function to demonstrate circular orbit calculations"""
    
    print("CIRCULAR ORBIT POSITION AND TIME CALCULATOR")
    print("=" * 50)
    
    # Create orbit objects
    leo_satellite = create_leo_satellite()
    moon = create_moon_orbit()
    
    # Generate orbit data for one LEO satellite orbit with 1-minute intervals
    print(f"\nGenerating orbit data for one LEO satellite orbit ({leo_satellite.period/60:.1f} minutes) with 1-minute intervals...")
    leo_data = leo_satellite.generate_orbit_data(time_interval_s=60, num_orbits=1)
    
    # For porkchop plots, we need much more Moon data to cover maximum time of flight
    # Simulate enough Moon orbit to cover departure duration + max time of flight + buffer
    max_time_hours = 24 * 5  # 5 days should cover max departure time + max TOF + buffer
    moon_orbits_needed = max_time_hours / (moon.period / 3600)  # Convert to number of Moon orbits
    print(f"Generating extended Moon orbit data for {max_time_hours} hours ({moon_orbits_needed:.3f} Moon orbits)...")
    moon_data = moon.generate_orbit_data(time_interval_s=60, num_orbits=moon_orbits_needed)
    
    print(f"LEO data points generated: {len(leo_data['time_min'])}")
    print(f"Moon data points generated: {len(moon_data['time_min'])}")
    
    # Display comparison table
    display_orbit_comparison(leo_data, moon_data)
    
    # Show specific positions at key times
    print(f"\n{'='*80}")
    print("KEY ORBITAL POSITIONS")
    print(f"{'='*80}")
    
    key_times = [0, 10, 20, 30, 44, 60, 80, 88]  # minutes - focus on one orbit duration
    
    print(f"\n{'Time':<15} {'LEO Position':<25} {'Moon Position':<25}")
    print("-" * 65)
    
    for t_min in key_times:
        t_sec = t_min * 60
        
        # LEO position
        leo_theta_rad, leo_theta_deg, leo_x, leo_y = leo_satellite.get_position_at_time(t_sec)
        # Moon position  
        moon_theta_rad, moon_theta_deg, moon_x, moon_y = moon.get_position_at_time(t_sec)
        
        if t_min < 60:
            time_str = f"{t_min} min"
        else:
            time_str = f"{t_min/60:.1f} hours"
            
        leo_pos_str = f"θ={leo_theta_deg%360:.1f}°"
        moon_pos_str = f"θ={moon_theta_deg%360:.3f}°"
        
        print(f"{time_str:<15} {leo_pos_str:<25} {moon_pos_str:<25}")
    
    # Plot the orbits - show FULL simulation time
    print(f"\nGenerating plots...")
    # Calculate maximum available time from both datasets
    max_leo_time = np.max(leo_data['time_h'])
    max_moon_time = np.max(moon_data['time_h'])
    full_sim_time = max(max_leo_time, max_moon_time)
    
    print(f"📊 Plotting full simulation time: {full_sim_time:.1f} hours ({full_sim_time/24:.1f} days)")
    plot_orbits(leo_data, moon_data, time_hours=full_sim_time)
    plot_angular_positions(leo_data, moon_data, time_hours=full_sim_time)
    
    # Save data to arrays for further analysis
    print(f"\nData arrays available for analysis:")
    print(f"  leo_data['time_min']: Time in minutes")
    print(f"  leo_data['theta_deg']: LEO angular position in degrees")
    print(f"  leo_data['x_km'], leo_data['y_km']: LEO Cartesian coordinates")
    print(f"  moon_data['theta_deg']: Moon angular position in degrees") 
    print(f"  moon_data['x_km'], moon_data['y_km']: Moon Cartesian coordinates")
    
    return leo_data, moon_data

if __name__ == "__main__":
    leo_data, moon_data = main()
