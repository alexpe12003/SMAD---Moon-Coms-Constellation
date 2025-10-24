#!/usr/bin/env python3
"""
Debug the return values from orbit conversion function
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from operations.trajectory_calculations import lunar_trajectory_calculations
from operations.lunar_operations import lunar_soi_calculations, hyperbolic_to_elliptical_conversion

def debug_return_values():
    """Check what the orbit conversion function actually returns"""
    print("DEBUGGING ORBIT CONVERSION RETURN VALUES")
    print("=" * 50)
    
    # Use option 3 parameters
    R0, V0, gamma0, lambda1 = 1.050, 1.372, 10.0, 229.1
    print(f"Testing with: R0={R0}, V0={V0}, γ0={gamma0}°, λ1={lambda1}°")
    
    try:
        geo_results = lunar_trajectory_calculations(R0, V0, gamma0, lambda1, verbose=False)
        lunar_results = lunar_soi_calculations(
            geo_results['r1'], geo_results['v1'], geo_results['phi1_deg'],
            lambda1, geo_results['gamma1_deg'], verbose=False
        )
        
        print(f"\nNatural perigee altitude: {lunar_results['hp']:.1f} km")
        
        # Get orbit conversion results
        orbit_results = hyperbolic_to_elliptical_conversion(lunar_results, verbose=False)
        
        print(f"\nOrbit conversion results:")
        print(f"=" * 30)
        
        # Print all keys and values
        for key, value in orbit_results.items():
            if isinstance(value, (int, float)):
                print(f"  {key}: {value:.3f}")
            else:
                print(f"  {key}: {value}")
        
        # Check specific fields
        print(f"\n🔍 KEY FIELDS:")
        print(f"  orbit_type: {orbit_results.get('orbit_type', 'NOT FOUND')}")
        print(f"  final_hp: {orbit_results.get('final_hp', 'NOT FOUND')}")
        print(f"  final_ha: {orbit_results.get('final_ha', 'NOT FOUND')}")
        print(f"  hp_circular (legacy): {orbit_results.get('hp_circular', 'NOT FOUND')}")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    debug_return_values()