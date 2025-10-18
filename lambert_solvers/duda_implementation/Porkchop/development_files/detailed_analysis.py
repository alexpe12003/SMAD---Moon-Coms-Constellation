"""
Detailed Circular Orbit Analysis
================================
This script provides detailed mathematical analysis of the circular orbits
defined in porkchop.py, including derivations and additional calculations.

Mathematical Background:
- For circular orbits: v = √(μ/r)
- Orbital period: T = 2π√(r³/μ)
- Angular velocity: ω = 2π/T = √(μ/r³)
- Position: θ(t) = ωt (starting from θ=0)
- Cartesian coordinates: x(t) = r·cos(θ(t)), y(t) = r·sin(θ(t))
"""

import numpy as np
import matplotlib.pyplot as plt
import math
from porkchop import CircularOrbit, create_leo_satellite, create_moon_orbit

def mathematical_analysis():
    """Provide detailed mathematical analysis of the orbits"""
    
    print("MATHEMATICAL ANALYSIS OF CIRCULAR ORBITS")
    print("=" * 50)
    
    # Constants
    MU_EARTH = 398600.4418  # km³/s²
    EARTH_RADIUS = 6378     # km
    MOON_DISTANCE = 384400  # km
    
    print(f"\nGiven Constants:")
    print(f"  μ_Earth = {MU_EARTH} km³/s²")
    print(f"  R_Earth = {EARTH_RADIUS} km")
    print(f"  Earth-Moon distance = {MOON_DISTANCE} km")
    
    # LEO calculations
    print(f"\nLEO SATELLITE (200km altitude):")
    print(f"{'='*40}")
    r_leo = EARTH_RADIUS + 200
    print(f"  Orbital radius: r = R_Earth + altitude = {EARTH_RADIUS} + 200 = {r_leo} km")
    
    v_leo = math.sqrt(MU_EARTH / r_leo)
    print(f"  Orbital velocity: v = √(μ/r) = √({MU_EARTH}/{r_leo}) = {v_leo:.3f} km/s")
    
    T_leo = 2 * math.pi * math.sqrt(r_leo**3 / MU_EARTH)
    print(f"  Orbital period: T = 2π√(r³/μ) = 2π√({r_leo}³/{MU_EARTH}) = {T_leo:.0f} s")
    print(f"                    = {T_leo/60:.1f} minutes = {T_leo/3600:.2f} hours")
    
    omega_leo = math.sqrt(MU_EARTH / r_leo**3)
    print(f"  Angular velocity: ω = √(μ/r³) = √({MU_EARTH}/{r_leo}³) = {omega_leo:.6f} rad/s")
    print(f"                      = {math.degrees(omega_leo):.6f} deg/s")
    
    # Moon calculations
    print(f"\nMOON ORBIT:")
    print(f"{'='*40}")
    r_moon = MOON_DISTANCE
    print(f"  Orbital radius: r = {r_moon} km")
    
    v_moon = math.sqrt(MU_EARTH / r_moon)
    print(f"  Orbital velocity: v = √(μ/r) = √({MU_EARTH}/{r_moon}) = {v_moon:.3f} km/s")
    
    T_moon = 2 * math.pi * math.sqrt(r_moon**3 / MU_EARTH)
    print(f"  Orbital period: T = 2π√(r³/μ) = 2π√({r_moon}³/{MU_EARTH}) = {T_moon:.0f} s")
    print(f"                    = {T_moon/60:.1f} minutes = {T_moon/3600:.1f} hours = {T_moon/86400:.1f} days")
    
    omega_moon = math.sqrt(MU_EARTH / r_moon**3)
    print(f"  Angular velocity: ω = √(μ/r³) = √({MU_EARTH}/{r_moon}³) = {omega_moon:.9f} rad/s")
    print(f"                      = {math.degrees(omega_moon):.6f} deg/s")
    
    # Comparison
    print(f"\nCOMPARATIVE ANALYSIS:")
    print(f"{'='*40}")
    print(f"  Radius ratio (Moon/LEO): {r_moon/r_leo:.1f}")
    print(f"  Velocity ratio (LEO/Moon): {v_leo/v_moon:.1f}")
    print(f"  Period ratio (Moon/LEO): {T_moon/T_leo:.1f}")
    print(f"  Angular velocity ratio (LEO/Moon): {omega_leo/omega_moon:.1f}")
    
    return {
        'leo': {'r': r_leo, 'v': v_leo, 'T': T_leo, 'omega': omega_leo},
        'moon': {'r': r_moon, 'v': v_moon, 'T': T_moon, 'omega': omega_moon}
    }

def position_equations_demo():
    """Demonstrate the position equations at specific times"""
    
    print(f"\nPOSITION EQUATIONS DEMONSTRATION")
    print(f"{'='*50}")
    
    # Create orbits
    leo = create_leo_satellite()
    moon = create_moon_orbit()
    
    print(f"\nPosition equations:")
    print(f"  θ(t) = ω × t")
    print(f"  x(t) = r × cos(θ(t))")
    print(f"  y(t) = r × sin(θ(t))")
    print(f"\nwhere t is time in seconds from initial position (θ=0)")
    
    # Demonstrate at specific times
    times_min = [0, 15, 30, 45, 60, 90, 120]
    
    print(f"\n{'Time':<8} {'LEO θ(°)':<10} {'LEO x(km)':<12} {'LEO y(km)':<12} {'Moon θ(°)':<12} {'Moon x(km)':<14} {'Moon y(km)':<14}")
    print("-" * 90)
    
    for t_min in times_min:
        t_sec = t_min * 60
        
        # LEO calculations
        theta_leo_rad = leo.angular_velocity * t_sec
        theta_leo_deg = math.degrees(theta_leo_rad)
        x_leo = leo.radius * math.cos(theta_leo_rad)
        y_leo = leo.radius * math.sin(theta_leo_rad)
        
        # Moon calculations
        theta_moon_rad = moon.angular_velocity * t_sec
        theta_moon_deg = math.degrees(theta_moon_rad)
        x_moon = moon.radius * math.cos(theta_moon_rad)
        y_moon = moon.radius * math.sin(theta_moon_rad)
        
        print(f"{t_min:<8} {theta_leo_deg%360:<10.1f} {x_leo:<12.0f} {y_leo:<12.0f} {theta_moon_deg:<12.3f} {x_moon:<14.0f} {y_moon:<14.0f}")

def synodic_period_analysis():
    """Calculate synodic period between LEO satellite and Moon"""
    
    print(f"\nSYNODIC PERIOD ANALYSIS")
    print(f"{'='*50}")
    
    leo = create_leo_satellite()
    moon = create_moon_orbit()
    
    # Synodic period calculation
    # 1/T_synodic = |1/T_leo - 1/T_moon|
    synodic_period = abs(1/(1/leo.period - 1/moon.period))
    
    print(f"\nSynodic period calculation:")
    print(f"  The synodic period is the time between successive alignments")
    print(f"  Formula: 1/T_synodic = |1/T_leo - 1/T_moon|")
    print(f"  1/T_synodic = |1/{leo.period:.0f} - 1/{moon.period:.0f}|")
    print(f"  1/T_synodic = |{1/leo.period:.6f} - {1/moon.period:.6f}|")
    print(f"  1/T_synodic = {abs(1/leo.period - 1/moon.period):.6f}")
    print(f"  T_synodic = {synodic_period:.0f} seconds = {synodic_period/3600:.2f} hours")
    
    # Number of LEO orbits in synodic period
    leo_orbits_in_synodic = synodic_period / leo.period
    moon_orbits_in_synodic = synodic_period / moon.period
    
    print(f"\nDuring one synodic period:")
    print(f"  LEO completes {leo_orbits_in_synodic:.1f} orbits")
    print(f"  Moon completes {moon_orbits_in_synodic:.3f} orbits")
    
    return synodic_period

def create_summary_table():
    """Create a comprehensive summary table"""
    
    print(f"\nCOMPREHENSIVE ORBITAL SUMMARY")
    print(f"{'='*70}")
    
    # Parameters
    params = {
        'Parameter': ['Radius (km)', 'Altitude (km)', 'Velocity (km/s)', 
                     'Period (hours)', 'Period (days)', 'Angular velocity (deg/s)',
                     'Circumference (km)', 'Angular velocity (deg/min)'],
        'LEO Satellite': [],
        'Moon': []
    }
    
    leo = create_leo_satellite()
    moon = create_moon_orbit()
    
    # LEO values
    params['LEO Satellite'] = [
        f"{leo.radius:.0f}",
        f"{leo.radius - 6378:.0f}",
        f"{leo.velocity:.3f}",
        f"{leo.period/3600:.2f}",
        f"{leo.period/86400:.3f}",
        f"{math.degrees(leo.angular_velocity):.6f}",
        f"{2*math.pi*leo.radius:.0f}",
        f"{math.degrees(leo.angular_velocity)*60:.3f}"
    ]
    
    # Moon values
    params['Moon'] = [
        f"{moon.radius:.0f}",
        f"{moon.radius - 6378:.0f}",
        f"{moon.velocity:.3f}",
        f"{moon.period/3600:.1f}",
        f"{moon.period/86400:.1f}",
        f"{math.degrees(moon.angular_velocity):.6f}",
        f"{2*math.pi*moon.radius:.0f}",
        f"{math.degrees(moon.angular_velocity)*60:.6f}"
    ]
    
    # Print table
    print(f"\n{'Parameter':<25} {'LEO Satellite':<20} {'Moon':<20}")
    print("-" * 65)
    for i, param in enumerate(params['Parameter']):
        print(f"{param:<25} {params['LEO Satellite'][i]:<20} {params['Moon'][i]:<20}")

def generate_position_arrays(duration_type='leo_orbit', interval_minutes=1):
    """Generate position arrays for both orbits"""
    
    leo = create_leo_satellite()
    
    if duration_type == 'leo_orbit':
        duration_hours = leo.period / 3600  # One LEO orbit in hours
        duration_desc = f"one LEO orbit ({leo.period/60:.1f} minutes)"
    else:
        duration_hours = 24  # 24 hours
        duration_desc = "24 hours"
    
    print(f"\nGENERATING POSITION ARRAYS")
    print(f"{'='*50}")
    print(f"Duration: {duration_desc}")
    print(f"Interval: {interval_minutes} minute(s)")
    
    moon = create_moon_orbit()
    
    # Time array
    times_min = np.arange(0, duration_hours * 60 + interval_minutes, interval_minutes)
    times_sec = times_min * 60
    
    # LEO arrays
    leo_theta_rad = leo.angular_velocity * times_sec
    leo_theta_deg = np.degrees(leo_theta_rad)
    leo_x = leo.radius * np.cos(leo_theta_rad)
    leo_y = leo.radius * np.sin(leo_theta_rad)
    
    # Moon arrays
    moon_theta_rad = moon.angular_velocity * times_sec
    moon_theta_deg = np.degrees(moon_theta_rad)
    moon_x = moon.radius * np.cos(moon_theta_rad)
    moon_y = moon.radius * np.sin(moon_theta_rad)
    
    print(f"Generated {len(times_min)} data points for each orbit")
    
    # Sample output
    print(f"\nSample data (first 5 points):")
    print(f"{'Time(min)':<10} {'LEO θ(°)':<10} {'Moon θ(°)':<12}")
    print("-" * 32)
    for i in range(5):
        print(f"{times_min[i]:<10.0f} {leo_theta_deg[i]%360:<10.1f} {moon_theta_deg[i]:<12.6f}")
    
    return {
        'time_min': times_min,
        'time_sec': times_sec,
        'leo': {
            'theta_rad': leo_theta_rad,
            'theta_deg': leo_theta_deg,
            'x_km': leo_x,
            'y_km': leo_y
        },
        'moon': {
            'theta_rad': moon_theta_rad,
            'theta_deg': moon_theta_deg,
            'x_km': moon_x,
            'y_km': moon_y
        }
    }

def main():
    """Main analysis function"""
    
    print("DETAILED CIRCULAR ORBIT ANALYSIS")
    print("=" * 50)
    
    # Mathematical analysis
    math_results = mathematical_analysis()
    
    # Position equations demonstration
    position_equations_demo()
    
    # Synodic period analysis
    synodic_period = synodic_period_analysis()
    
    # Summary table
    create_summary_table()
    
    # Generate position arrays
    data = generate_position_arrays(duration_type='leo_orbit', interval_minutes=1)
    
    print(f"\nANALYSIS COMPLETE")
    print(f"{'='*50}")
    print(f"All orbital parameters calculated and arrays generated.")
    print(f"Data available for further analysis and plotting.")
    
    return data

if __name__ == "__main__":
    data = main()