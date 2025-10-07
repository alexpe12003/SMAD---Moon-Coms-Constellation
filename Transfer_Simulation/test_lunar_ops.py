#!/usr/bin/env python3
"""
Test script to validate lunar_operations.py against Transfer.py
"""

import sys
import os

# Add the current directory to Python path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from lunar_operations import lunar_soi_calculations, hyperbolic_to_elliptical_conversion

def test_lunar_operations():
    """Test the lunar operations with known values"""
    
    print("TESTING LUNAR OPERATIONS MODULE")
    print("=" * 50)
    
    # Test parameters (from successful Transfer.py run)
    r1 = 51.5305  # DU
    v1 = 0.1312   # DU/TU  
    phi1_deg = 77.70  # degrees
    lambda1_deg = 30.0  # degrees
    gamma1_deg = 5.79   # degrees
    
    print(f"Test inputs:")
    print(f"  r1 = {r1} DU")
    print(f"  v1 = {v1} DU/TU")
    print(f"  phi1 = {phi1_deg}°")
    print(f"  lambda1 = {lambda1_deg}°")
    print(f"  gamma1 = {gamma1_deg}°")
    print()
    
    # Test lunar SOI calculations
    print("1. Testing lunar_soi_calculations()...")
    lunar_results = lunar_soi_calculations(r1, v1, phi1_deg, lambda1_deg, gamma1_deg, verbose=True)
    
    print("\n" + "=" * 50)
    print("2. Testing hyperbolic_to_elliptical_conversion()...")
    elliptical_results = hyperbolic_to_elliptical_conversion(lunar_results, target_perigee_altitude_km=600, verbose=True)
    
    print("\n" + "=" * 50)
    print("SUMMARY OF RESULTS:")
    print(f"  Lunar periapsis altitude: {lunar_results['hp']:.0f} km")
    print(f"  Total delta-V required: {elliptical_results['total_delta_v']:.3f} km/s")
    print(f"  Hyperbolic to elliptical: {abs(elliptical_results['delta_v_conversion']):.3f} km/s")
    print(f"  Elliptical to circular: {abs(elliptical_results['delta_v_circularization']):.3f} km/s")

if __name__ == "__main__":
    test_lunar_operations()