"""
User interface and input/output functions for lunar transfer analysis
"""

from config import *


def get_user_input():
    """
    Get mission parameters including lambda1
    Returns a dictionary with the input values
    """
    print("Please enter the following parameters:")
    print("-" * 50)
    
    variables = {}
    
    # Set predefined mission parameters
    variables['R0'] = DEFAULT_PARKING_ORBIT_RADIUS  # DU - Both parking orbit radius AND transfer trajectory perigee
    variables['V0'] = DEFAULT_TRANSFER_VELOCITY  # DU/TU - Transfer trajectory velocity at perigee
    variables['gamma0'] = DEFAULT_FLIGHT_PATH_ANGLE  # degrees - Flight path angle at perigee
    
    print(f"Using predefined mission parameters:")
    print(f"Parking orbit radius = Transfer perigee (R0) = {variables['R0']} DU")
    print(f"Parking orbit altitude = {(variables['R0']-1)*DU_TO_KM:.0f} km")
    print(f"Transfer trajectory velocity (V0) = {variables['V0']} DU/TU")
    print(f"Transfer flight path angle (gamma0) = {variables['gamma0']}°")
    print()
    
    # Only ask for lambda1 from user
    while True:
        try:
            value = float(input("Enter lambda1 (Point where geocentric trajectory crosses lunar sphere of influence in degrees): "))
            variables['lambda1'] = value
            break
        except ValueError:
            print("Error: Please enter a valid number.")
            continue
    
    return variables


def display_organized_mission_summary(inputs, departure_results, geo_results, lunar_results, soi_transit_results, elliptical_results=None):
    """
    Display a comprehensive, organized summary of the complete lunar transfer mission
    
    Parameters:
    - inputs: User input parameters
    - departure_results: Earth departure analysis results
    - geo_results: Geocentric trajectory results
    - lunar_results: Lunar SOI trajectory results
    - soi_transit_results: SOI transit time results
    - elliptical_results: Orbit conversion results (optional)
    """
    print("\n" + "=" * 80)
    print("COMPLETE LUNAR TRANSFER MISSION ANALYSIS")
    print("=" * 80)
    
    # Mission Parameters
    print(f"\n📋 MISSION PARAMETERS:")
    print(f"   Parking Orbit = Transfer Perigee: {departure_results['parking_altitude_km']:.0f} km altitude")
    print(f"   Transfer Velocity (V0): {inputs['V0']} DU/TU")
    print(f"   Flight Path Angle (γ0): {inputs['gamma0']}°")
    print(f"   Lunar SOI Crossing (λ1): {inputs['lambda1']}°")
    
    # Step 0: Earth Departure
    print(f"\n🚀 STEP 0: EARTH DEPARTURE")
    print(f"   Parking Orbit Velocity: {departure_results['V_circular_parking_kms']:.3f} km/s")
    print(f"   Transfer Velocity at Departure: {departure_results['V0_kms']:.3f} km/s")
    print(f"   Departure Delta-V: {departure_results['delta_v_departure_kms']:.3f} km/s")
    print(f"   Maneuver Location: {departure_results['parking_altitude_km']:.0f} km altitude")
    
    # Step 1: Earth Departure Trajectory
    print(f"\n🛰️ STEP 1: EARTH DEPARTURE TRAJECTORY")
    trajectory_type = "Elliptical" if departure_results['specific_energy'] < 0 else "Hyperbolic"
    print(f"   Trajectory Type: {trajectory_type} (e = {geo_results['e']:.3f})")
    print(f"   Specific Energy: {departure_results['specific_energy']:.4f} DU²/TU²")
    print(f"   Angular Momentum: {geo_results['h']:.3f} DU²/TU")
    print(f"   Semi-major Axis: {geo_results['a']:.2f} DU")
    
    # Step 2: Earth-Moon Transit
    print(f"\n🌍→🌙 STEP 2: EARTH-MOON TRANSIT")
    print(f"   Time of Flight: {geo_results['tof_hours']:.1f} hours ({geo_results['tof_TU']:.2f} TU)")
    print(f"   Distance at Moon's SOI: {geo_results['r1']:.2f} DU")
    print(f"   Velocity at Moon's SOI: {geo_results['v1']:.4f} DU/TU")
    
    # Step 3: Lunar SOI Approach
    print(f"\n🌙 STEP 3: LUNAR SPHERE OF INFLUENCE APPROACH")
    print(f"   Relative Velocity to Moon: {lunar_results['v2_kms']:.3f} km/s")
    print(f"   Flight Path Angle: {lunar_results['epsilon2_deg']:.1f}°")
    print(f"   Trajectory Type: Hyperbolic (e = {lunar_results['e_lunar']:.3f})")
    print(f"   SOI Entry to Perigee Time: {soi_transit_results['soi_transit_time_hours']:.2f} hours")
    print(f"   Natural Flyby Altitude: {lunar_results['hp']:.0f} km above surface")
    
    # Step 4: Lunar Orbit Insertion (if elliptical conversion provided)
    if elliptical_results:
        print(f"\n🛰️ STEP 4: LUNAR ORBIT INSERTION")
        print(f"   Maneuver Location: {elliptical_results['hp_hyperbolic']:.0f} km altitude")
        print(f"   Delta-V Required: {elliptical_results['delta_v_conversion']:.3f} km/s")
        print(f"   Resulting Orbit: {elliptical_results['hp_hyperbolic']:.0f} × {elliptical_results['target_perigee_altitude_km']:.0f} km")
        print(f"   Orbital Period: {elliptical_results['orbital_period_hours']:.2f} hours")
        
        # Step 5: Orbit Circularization
        print(f"\n🔄 STEP 5: ORBIT CIRCULARIZATION (Optional)")
        print(f"   Circularization at: {elliptical_results['target_perigee_altitude_km']:.0f} km altitude")
        print(f"   Delta-V Required: {elliptical_results['delta_v_circularization']:.3f} km/s")
    
    # Mission Summary
    print(f"\n📊 COMPLETE MISSION SUMMARY:")
    print(f"   Earth Departure Delta-V: {departure_results['delta_v_departure_kms']:.3f} km/s")
    if elliptical_results:
        print(f"   Lunar Operations Delta-V: {elliptical_results['total_delta_v']:.3f} km/s")
        total_mission_dv = departure_results['delta_v_departure_kms'] + elliptical_results['total_delta_v']
        print(f"   TOTAL MISSION DELTA-V: {total_mission_dv:.3f} km/s")
    
    print(f"   Earth-Moon Transit Time: {geo_results['tof_hours']:.1f} hours")
    print(f"   SOI Transit Time: {soi_transit_results['soi_transit_time_hours']:.2f} hours")
    total_mission_time = geo_results['tof_hours'] + soi_transit_results['soi_transit_time_hours']
    print(f"   TOTAL MISSION DURATION: {total_mission_time:.1f} hours")
    
    if elliptical_results:
        print(f"   Final Orbit: {elliptical_results['target_perigee_altitude_km']:.0f} km circular")
    
    print("=" * 80)


def display_analysis_menu():
    """
    Display the main analysis menu and get user choice
    
    Returns:
    - String: User's choice ('1' or '2')
    """
    print("LUNAR TRANSFER TRAJECTORY ANALYSIS PROGRAM")
    print("=" * 60)
    
    # Ask user what analysis to perform
    print("Choose analysis type:")
    print("1. Single calculation with user-defined lambda1")
    print("2. Parametric study: lambda1 from 0-360° (5° steps)")
    
    while True:
        try:
            choice = input("Enter choice (1 or 2): ").strip()
            if choice in ['1', '2']:
                return choice
            else:
                print("Please enter 1 or 2")
        except:
            print("Please enter 1 or 2")