#!/usr/bin/env python3
"""
Quick test of option 3 with corrected logic
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from operations.trajectory_calculations import lunar_trajectory_calculations
from operations.lunar_operations import lunar_soi_calculations, hyperbolic_to_elliptical_conversion

def test_option3():
    """Test option 3 parameters with corrected logic"""
    print("TESTING OPTION 3 PARAMETERS")
    print("=" * 40)
    
    # Current option 3 parameters
    R0, V0, gamma0, lambda1 = 1.050, 1.372, 10.0, 229.1
    print(f"Parameters: R0={R0}, V0={V0}, γ0={gamma0}°, λ1={lambda1}°")
    
    try:
        geo_results = lunar_trajectory_calculations(R0, V0, gamma0, lambda1, verbose=False)
        lunar_results = lunar_soi_calculations(
            geo_results['r1'], geo_results['v1'], geo_results['phi1_deg'],
            lambda1, geo_results['gamma1_deg'], verbose=False
        )
        
        perigee_alt = lunar_results['hp']
        print(f"\nNatural perigee altitude: {perigee_alt:.1f} km")
        print(f"Target apoapsis altitude: 1500 km")
        
        # Test orbit conversion
        orbit_results = hyperbolic_to_elliptical_conversion(lunar_results, verbose=True)
        
        if not orbit_results.get('error'):
            print(f"\n🎯 FINAL RESULT:")
            print(f"   Orbit type: {orbit_results['orbit_type']}")
            print(f"   Total Delta-V: {orbit_results['total_delta_v']:.3f} km/s")
            
            if orbit_results['orbit_type'] == 'elliptical':
                print(f"   Perigee altitude: {orbit_results['final_hp']:.1f} km")  
                print(f"   Apoapsis altitude: {orbit_results['final_ha']:.1f} km")
                
                # Verify correctness
                if perigee_alt < 1500 and orbit_results['final_ha'] == 1500:
                    print(f"   ✅ CORRECT: perigee < 1500km, apoapsis = 1500km")
                else:
                    print(f"   ❌ LOGIC ERROR")
                    
        else:
            print(f"❌ Error: {orbit_results['error_message']}")
            
    except Exception as e:
        print(f"❌ Exception: {e}")

if __name__ == "__main__":
    test_option3()