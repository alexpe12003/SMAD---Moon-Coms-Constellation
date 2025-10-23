#!/usr/bin/env python3
"""
Test the corrected 1500km apoapsis logic
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from operations.trajectory_calculations import lunar_trajectory_calculations
from operations.lunar_operations import lunar_soi_calculations, hyperbolic_to_elliptical_conversion

def test_corrected_logic():
    """Test the corrected 1500km apoapsis logic"""
    print("TESTING CORRECTED 1500KM APOAPSIS LOGIC")
    print("=" * 60)
    
    test_cases = [
        # Test case 1: Known safe parameters (perigee ~1084 km < 1500 km)
        {"name": "Safe trajectory (perigee < 1500km)", "params": (1.05, 1.372, 10.0, 229.1)},
        
        # Test case 2: Try to find a case with perigee > 1500km  
        {"name": "High perigee trajectory", "params": (1.05, 1.372, 10.0, 180.0)},
        
        # Test case 3: Collision parameters (should be rejected)
        {"name": "Collision trajectory", "params": (1.05, 1.372, 29.5, 36.4)},
    ]
    
    for i, test_case in enumerate(test_cases, 1):
        print(f"\nTest {i}: {test_case['name']}")
        print("-" * 40)
        
        R0, V0, gamma0, lambda1 = test_case['params']
        print(f"Parameters: R0={R0}, V0={V0}, γ0={gamma0}°, λ1={lambda1}°")
        
        try:
            # Geocentric trajectory
            geo_results = lunar_trajectory_calculations(R0, V0, gamma0, lambda1, verbose=False)
            
            # Lunar SOI calculations  
            lunar_results = lunar_soi_calculations(
                geo_results['r1'], geo_results['v1'], geo_results['phi1_deg'],
                lambda1, geo_results['gamma1_deg'], verbose=False
            )
            
            perigee_alt = lunar_results['hp']
            print(f"Natural perigee altitude: {perigee_alt:.1f} km")
            
            # Orbit conversion with new logic
            orbit_results = hyperbolic_to_elliptical_conversion(lunar_results, verbose=False)
            
            if orbit_results.get('error') == 'COLLISION_TRAJECTORY':
                print(f"✅ Collision correctly detected and rejected")
                
            elif orbit_results.get('orbit_type') == 'elliptical':
                print(f"✅ Elliptical orbit created:")
                print(f"   Perigee: {orbit_results['final_hp']:.1f} km")
                print(f"   Apoapsis: {orbit_results['final_ha']:.1f} km")
                print(f"   Delta-V: {orbit_results['total_delta_v']:.3f} km/s")
                
                # Verify logic
                if perigee_alt < 1500:
                    if orbit_results['final_ha'] == 1500:
                        print(f"✅ Correct: perigee < 1500km → apoapsis = 1500km")
                    else:
                        print(f"❌ Error: Expected apoapsis = 1500km, got {orbit_results['final_ha']:.1f}km")
                
            elif orbit_results.get('orbit_type') == 'circular':
                print(f"✅ Circular orbit created:")
                print(f"   Altitude: {orbit_results['final_hp']:.1f} km")
                print(f"   Delta-V: {orbit_results['total_delta_v']:.3f} km/s")
                
                # Verify logic
                if perigee_alt >= 1500:
                    print(f"✅ Correct: perigee ≥ 1500km → circularized")
                else:
                    print(f"❌ Error: Expected elliptical orbit for perigee < 1500km")
                    
            else:
                print(f"❌ Unknown orbit type: {orbit_results.get('orbit_type')}")
                
        except Exception as e:
            print(f"❌ Error: {e}")
    
    print(f"\n" + "=" * 60)
    print("SUMMARY:")
    print("✅ Logic should be: perigee < 1500km → ellipse with apoapsis at 1500km")
    print("✅ Logic should be: perigee ≥ 1500km → circularize at natural perigee") 
    print("✅ Collision trajectories should be rejected")

if __name__ == "__main__":
    test_corrected_logic()