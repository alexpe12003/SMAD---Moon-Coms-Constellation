"""
Position and Time Data Access Script
===================================
This script demonstrates how to access the position and time data for:
1. LEO satellite in 200km circular orbit
2. Moon in circular orbit at average Eart    print(f"\n" + "=" * 50)
    print("DATA GENERATION COMPLETE")
    print("=" * 50)
    print("You now have:")
    print("1. leo_data - dictionary with LEO satellite position/time data")
    print("2. moon_data - dictionary with Moon position/time data")
    print("3. All data saved to 'orbit_data.txt'")
    print("\nBoth orbits are coplanar, 2D, with theta=0 at initial position")
    print("Time intervals: 1 minute")
    print(f"Duration: {leo_satellite.period/60:.1f} minutes (1 complete LEO orbit)") print("Both orbits are coplanar, 2D, with theta=0 at initial position")istance

Both orbits are coplanar and 2D, with position defined by angle θ
where θ=0 at initial position, with 1-minute intervals.
"""

import numpy as np
import math
from porkchop import create_leo_satellite, create_moon_orbit

def generate_orbit_data_tables():
    """Generate and display orbit data in tabular format"""
    
    print("ORBIT DATA GENERATION - 1 MINUTE INTERVALS")
    print("=" * 60)
    
    # Create orbit objects
    leo_satellite = create_leo_satellite()
    moon = create_moon_orbit()
    
    # Time parameters
    time_interval_minutes = 1
    duration_minutes = leo_satellite.period / 60  # One LEO satellite orbit
    
    # Generate time array (in minutes)
    time_minutes = np.arange(0, duration_minutes + time_interval_minutes, time_interval_minutes)
    time_seconds = time_minutes * 60
    
    print(f"\nData generation parameters:")
    print(f"  Time interval: {time_interval_minutes} minute(s)")
    print(f"  Duration: {duration_minutes:.1f} minutes ({duration_minutes/60:.1f} hours)")
    print(f"  Total data points: {len(time_minutes)}")
    
    # LEO Satellite Data
    print(f"\nLEO SATELLITE DATA:")
    print(f"{'='*50}")
    
    leo_data = {
        'time_min': time_minutes,
        'time_sec': time_seconds,
        'theta_rad': [],
        'theta_deg': [],
        'x_km': [],
        'y_km': []
    }
    
    for t_sec in time_seconds:
        theta_rad, theta_deg, x_km, y_km = leo_satellite.get_position_at_time(t_sec)
        leo_data['theta_rad'].append(theta_rad)
        leo_data['theta_deg'].append(theta_deg)
        leo_data['x_km'].append(x_km)
        leo_data['y_km'].append(y_km)
    
    # Convert to numpy arrays
    for key in ['theta_rad', 'theta_deg', 'x_km', 'y_km']:
        leo_data[key] = np.array(leo_data[key])
    
    # Moon Data
    print(f"\nMOON DATA:")
    print(f"{'='*50}")
    
    moon_data = {
        'time_min': time_minutes,
        'time_sec': time_seconds,
        'theta_rad': [],
        'theta_deg': [],
        'x_km': [],
        'y_km': []
    }
    
    for t_sec in time_seconds:
        theta_rad, theta_deg, x_km, y_km = moon.get_position_at_time(t_sec)
        moon_data['theta_rad'].append(theta_rad)
        moon_data['theta_deg'].append(theta_deg)
        moon_data['x_km'].append(x_km)
        moon_data['y_km'].append(y_km)
    
    # Convert to numpy arrays
    for key in ['theta_rad', 'theta_deg', 'x_km', 'y_km']:
        moon_data[key] = np.array(moon_data[key])
    
    return leo_data, moon_data, leo_satellite

def display_data_samples(leo_data, moon_data, num_samples=20):
    """Display sample data points"""
    
    print(f"\nSAMPLE DATA POINTS (First {num_samples} points)")
    print("=" * 100)
    
    print(f"{'Time':<6} {'LEO θ':<8} {'LEO X':<10} {'LEO Y':<10} {'Moon θ':<10} {'Moon X':<12} {'Moon Y':<12}")
    print(f"{'(min)':<6} {'(deg)':<8} {'(km)':<10} {'(km)':<10} {'(deg)':<10} {'(km)':<12} {'(km)':<12}")
    print("-" * 100)
    
    for i in range(min(num_samples, len(leo_data['time_min']))):
        time_min = leo_data['time_min'][i]
        leo_theta = leo_data['theta_deg'][i] % 360  # Normalize to 0-360°
        leo_x = leo_data['x_km'][i]
        leo_y = leo_data['y_km'][i]
        moon_theta = moon_data['theta_deg'][i] % 360  # Normalize to 0-360°
        moon_x = moon_data['x_km'][i]
        moon_y = moon_data['y_km'][i]
        
        print(f"{time_min:<6.0f} {leo_theta:<8.1f} {leo_x:<10.0f} {leo_y:<10.0f} {moon_theta:<10.3f} {moon_x:<12.0f} {moon_y:<12.0f}")

def analyze_orbital_motion(leo_data, moon_data, leo_satellite):
    """Analyze key aspects of orbital motion"""
    
    print(f"\nORBITAL MOTION ANALYSIS")
    print("=" * 50)
    
    # Calculate angular velocities from data
    dt = (leo_data['time_sec'][1] - leo_data['time_sec'][0])  # Time step in seconds
    
    leo_angular_vel = (leo_data['theta_rad'][1] - leo_data['theta_rad'][0]) / dt  # rad/s
    moon_angular_vel = (moon_data['theta_rad'][1] - moon_data['theta_rad'][0]) / dt  # rad/s
    
    print(f"\nAngular velocities (calculated from data):")
    print(f"  LEO satellite: {leo_angular_vel:.6f} rad/s = {math.degrees(leo_angular_vel):.3f} deg/s")
    print(f"  Moon: {moon_angular_vel:.9f} rad/s = {math.degrees(moon_angular_vel):.6f} deg/s")
    
    # Find when LEO completes orbits
    leo_complete_orbits = []
    for i, theta in enumerate(leo_data['theta_deg']):
        if i > 0 and theta % 360 < leo_data['theta_deg'][i-1] % 360:
            leo_complete_orbits.append(i)
    
    print(f"\nLEO satellite completes orbits at:")
    for i, idx in enumerate(leo_complete_orbits[:5]):  # Show first 5
        time_min = leo_data['time_min'][idx]
        print(f"  Orbit {i+1}: {time_min:.1f} minutes ({time_min/60:.2f} hours)")
    
    # Moon angular displacement in one LEO orbit
    moon_angular_displacement = moon_data['theta_deg'][-1]
    print(f"\nMoon angular displacement in one LEO orbit ({leo_satellite.period/60:.1f} minutes):")
    print(f"  {moon_angular_displacement:.3f} degrees")
    print(f"  {moon_angular_displacement/360:.6f} complete orbits")

def create_data_access_examples(leo_data, moon_data):
    """Show examples of how to access specific data"""
    
    print(f"\nDATA ACCESS EXAMPLES")
    print("=" * 50)
    
    print(f"\nHow to access data arrays:")
    print(f"  leo_data['time_min']   - Time in minutes")
    print(f"  leo_data['theta_deg']  - LEO angle in degrees") 
    print(f"  leo_data['x_km']       - LEO X coordinate in km")
    print(f"  leo_data['y_km']       - LEO Y coordinate in km")
    print(f"  moon_data['theta_deg'] - Moon angle in degrees")
    print(f"  moon_data['x_km']      - Moon X coordinate in km")
    print(f"  moon_data['y_km']      - Moon Y coordinate in km")
    
    print(f"\nExample: Position at t = 60 minutes:")
    idx_60min = 60  # 60th data point (60 minutes)
    if idx_60min < len(leo_data['time_min']):
        print(f"  Time: {leo_data['time_min'][idx_60min]} minutes")
        print(f"  LEO: theta = {leo_data['theta_deg'][idx_60min]:.1f}°, x = {leo_data['x_km'][idx_60min]:.0f} km, y = {leo_data['y_km'][idx_60min]:.0f} km")
        print(f"  Moon: theta = {moon_data['theta_deg'][idx_60min]:.3f}°, x = {moon_data['x_km'][idx_60min]:.0f} km, y = {moon_data['y_km'][idx_60min]:.0f} km")
    
    print(f"\nExample: Find when LEO satellite reaches theta = 90°:")
    for i, theta in enumerate(leo_data['theta_deg']):
        if abs((theta % 360) - 90) < 1:  # Within 1 degree of 90°
            print(f"  At t = {leo_data['time_min'][i]:.1f} min: theta = {theta%360:.1f}°")
            break
    
    print(f"\nExample: Calculate distance between LEO and Moon at t = 120 min:")
    idx_120min = 120
    if idx_120min < len(leo_data['time_min']):
        leo_x = leo_data['x_km'][idx_120min]
        leo_y = leo_data['y_km'][idx_120min]
        moon_x = moon_data['x_km'][idx_120min]
        moon_y = moon_data['y_km'][idx_120min]
        
        distance = math.sqrt((moon_x - leo_x)**2 + (moon_y - leo_y)**2)
        print(f"  LEO position: ({leo_x:.0f}, {leo_y:.0f}) km")
        print(f"  Moon position: ({moon_x:.0f}, {moon_y:.0f}) km")
        print(f"  Distance: {distance:.0f} km")

def save_data_to_file(leo_data, moon_data, filename="orbit_data.txt"):
    """Save data to a text file"""
    
    print(f"\nSAVING DATA TO FILE: {filename}")
    print("=" * 50)
    
    with open(filename, 'w', encoding='utf-8') as f:
        f.write("Circular Orbit Position and Time Data\n")
        f.write("=" * 50 + "\n")
        f.write("LEO Satellite: 200km altitude circular orbit\n")
        f.write("Moon: Circular orbit at average Earth-Moon distance\n")
        f.write("Time intervals: 1 minute\n")
        f.write("Position defined by angle theta, with theta=0 at initial position\n\n")
        
        f.write(f"{'Time(min)':<10} {'LEO_theta(deg)':<12} {'LEO_X(km)':<12} {'LEO_Y(km)':<12} {'Moon_theta(deg)':<14} {'Moon_X(km)':<14} {'Moon_Y(km)':<14}\n")
        f.write("-" * 100 + "\n")
        
        for i in range(len(leo_data['time_min'])):
            time_min = leo_data['time_min'][i]
            leo_theta = leo_data['theta_deg'][i] % 360
            leo_x = leo_data['x_km'][i]
            leo_y = leo_data['y_km'][i]
            moon_theta = moon_data['theta_deg'][i] % 360
            moon_x = moon_data['x_km'][i]
            moon_y = moon_data['y_km'][i]
            
            f.write(f"{time_min:<10.1f} {leo_theta:<12.3f} {leo_x:<12.1f} {leo_y:<12.1f} {moon_theta:<14.6f} {moon_x:<14.1f} {moon_y:<14.1f}\n")
    
    print(f"Data saved to {filename}")
    print(f"Total data points: {len(leo_data['time_min'])}")

def main():
    """Main function to generate and analyze orbit data"""
    
    print("POSITION AND TIME DATA ACCESS")
    print("=" * 50)
    
    # Generate orbit data
    leo_data, moon_data, leo_satellite = generate_orbit_data_tables()
    
    # Display sample data
    display_data_samples(leo_data, moon_data, num_samples=20)
    
    # Analyze orbital motion
    analyze_orbital_motion(leo_data, moon_data, leo_satellite)
    
    # Show data access examples
    create_data_access_examples(leo_data, moon_data)
    
    # Save data to file
    save_data_to_file(leo_data, moon_data)
    
    print(f"\n" + "=" * 50)
    print("DATA GENERATION COMPLETE")
    print("=" * 50)
    print("You now have:")
    print("1. leo_data - dictionary with LEO satellite position/time data")
    print("2. moon_data - dictionary with Moon position/time data")
    print("3. All data saved to 'orbit_data.txt'")
    print("\nBoth orbits are coplanar, 2D, with θ=0 at initial position")
    print("Time intervals: 1 minute")
    print(f"Duration: {leo_satellite.period/60:.1f} minutes (1 complete LEO orbit)")
    
    return leo_data, moon_data, leo_satellite

if __name__ == "__main__":
    leo_data, moon_data, leo_satellite = main()