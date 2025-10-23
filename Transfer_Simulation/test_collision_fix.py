#!/usr/bin/env python3
"""
Test the fixed collision detection
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

def test_collision_detection():
    """Test both collision and safe trajectories"""
    print("TESTING COLLISION DETECTION")
    print("=" * 50)
    
    from operations.trajectory_calculations import lunar_trajectory_calculations
    from operations.lunar_operations import lunar_soi_calculations, hyperbolic_to_elliptical_conversion
    
    # Test 1: Collision parameters (from original "optimal")
    print("Test 1: Collision parameters (γ0=29.5°, λ1=36.4°)")
    try:
        geo_results = lunar_trajectory_calculations(1.05, 1.372, 29.5, 36.4, verbose=False)
        lunar_results = lunar_soi_calculations(
            geo_results['r1'], geo_results['v1'], geo_results['phi1_deg'],
            36.4, geo_results['gamma1_deg'], verbose=False
        )
        print(f"  Perigee altitude: {lunar_results['hp']:.1f} km")
        
        elliptical_results = hyperbolic_to_elliptical_conversion(lunar_results, verbose=False)
        if elliptical_results.get('error') == 'COLLISION_TRAJECTORY':
            print("  ✅ Collision detected correctly!")
        else:
            print("  ❌ Collision not detected!")
            
    except Exception as e:
        print(f"  Error: {e}")
    
    # Test 2: Safe parameters
    print(f"\nTest 2: Safe parameters (γ0=10°, λ1=229.1°)")
    try:
        geo_results = lunar_trajectory_calculations(1.05, 1.372, 10.0, 229.1, verbose=False)
        lunar_results = lunar_soi_calculations(
            geo_results['r1'], geo_results['v1'], geo_results['phi1_deg'],
            229.1, geo_results['gamma1_deg'], verbose=False
        )
        print(f"  Perigee altitude: {lunar_results['hp']:.1f} km")
        
        elliptical_results = hyperbolic_to_elliptical_conversion(lunar_results, verbose=False)
        if elliptical_results.get('error') == 'COLLISION_TRAJECTORY':
            print("  ❌ False collision detection!")
        else:
            print(f"  ✅ Safe trajectory - Delta-V: {elliptical_results['total_delta_v']:.3f} km/s")
            
    except Exception as e:
        print(f"  Error: {e}")

if __name__ == "__main__":
    test_collision_detection()