"""
Earth-related calculations for lunar transfer missions
"""

import math
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import *


def calculate_earth_departure_delta_v(R0, V0, verbose=True):
    """
    Calculate the delta-V required for Earth departure from circular parking orbit
    to the transfer trajectory.
    
    Parameters:
    - R0: Parking orbit radius (also transfer perigee) in DU
    - V0: Transfer trajectory velocity at perigee in DU/TU
    - verbose: Whether to print detailed output
    
    Returns:
    - Dictionary containing departure analysis results
    """
    
    # Calculate circular parking orbit velocity
    V_circular_parking = math.sqrt(MU_EARTH / R0)  # DU/TU
    
    # Calculate departure delta-V
    delta_V_departure = V0 - V_circular_parking  # DU/TU
    
    # Convert to km/s
    V_circular_parking_kms = V_circular_parking * DU_TU_TO_KM_S
    V0_kms = V0 * DU_TU_TO_KM_S
    delta_V_departure_kms = delta_V_departure * DU_TU_TO_KM_S
    
    # Calculate parking orbit altitude
    parking_altitude_km = (R0 - 1) * DU_TO_KM
    
    # Calculate transfer trajectory properties
    # Using vis-viva equation: V² = μ(2/r - 1/a)
    # Rearranging for semi-major axis: a = 1/(2/r - V²/μ)
    specific_energy = V0**2 / 2 - MU_EARTH / R0  # DU²/TU²
    semi_major_axis = -MU_EARTH / (2 * specific_energy)  # DU
    
    if verbose:
        print("=" * 60)
        print("EARTH DEPARTURE ANALYSIS")
        print("=" * 60)
        print(f"Parking Orbit Analysis:")
        print(f"  Radius: {R0:.4f} DU ({parking_altitude_km:.0f} km altitude)")
        print(f"  Circular velocity: {V_circular_parking:.4f} DU/TU = {V_circular_parking_kms:.3f} km/s")
        print()
        print(f"Transfer Trajectory Analysis:")
        print(f"  Perigee radius: {R0:.4f} DU (same as parking orbit)")
        print(f"  Semi-major axis: {semi_major_axis:.2f} DU")
        print(f"  Trajectory type: {'Elliptical' if specific_energy < 0 else 'Hyperbolic'}")
        print(f"  Specific energy: {specific_energy:.4f} DU²/TU²")
        print(f"  Velocity at departure: {V0:.4f} DU/TU = {V0_kms:.3f} km/s")
        print()
        print(f"Departure Maneuver:")
        print(f"  Location: {parking_altitude_km:.0f} km altitude")
        print(f"  Delta-V required: {delta_V_departure:.4f} DU/TU = {delta_V_departure_kms:.3f} km/s")
        print(f"  Maneuver type: {'Prograde' if delta_V_departure > 0 else 'Retrograde'} burn ({'acceleration' if delta_V_departure > 0 else 'deceleration'})")
    
    return {
        'R0': R0,
        'V_circular_parking': V_circular_parking,
        'V_circular_parking_kms': V_circular_parking_kms,
        'V0': V0,
        'V0_kms': V0_kms,
        'delta_v_departure': delta_V_departure,
        'delta_v_departure_kms': delta_V_departure_kms,
        'parking_altitude_km': parking_altitude_km,
        'specific_energy': specific_energy,
        'semi_major_axis': semi_major_axis,
        'trajectory_type': 'Elliptical' if specific_energy < 0 else 'Hyperbolic'
    }