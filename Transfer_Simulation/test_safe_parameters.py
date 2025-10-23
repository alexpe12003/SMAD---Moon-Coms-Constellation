#!/usr/bin/env python3
"""
Test the updated optimal parameters to ensure they don't cause collision
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from operations.earth_operations import calculate_earth_departure_delta_v
from operations.trajectory_calculations import lunar_trajectory_calculations
from operations.lunar_operations import lunar_soi_calculations, calculate_lunar_soi_transit_time, hyperbolic_to_elliptical_conversion

def test_safe_parameters():
    """
    Test the safe parameters: R0=1.05, V0=1.372, γ0=0.0°, λ1=240.0°
    """
    print("Testing safe optimal parameters...")
    print("=" * 50)
    
    # Safe parameters
    R0 = 1.050
    V0 = 1.372
    gamma0 = 0.0
    lambda1 = 240.0
    
    print(f"Parameters: R0={R0}, V0={V0}, γ0={gamma0}°, λ1={lambda1}°")
    
    try:
        # Calculate complete mission
        departure_results = calculate_earth_departure_delta_v(R0, V0, verbose=False)
        geo_results = lunar_trajectory_calculations(R0, V0, gamma0, lambda1, verbose=False)
        lunar_results = lunar_soi_calculations(
            geo_results['r1'], 
            geo_results['v1'], 
            geo_results['phi1_deg'],
            lambda1,
            geo_results['gamma1_deg'],
            verbose=False
        )
        soi_transit_results = calculate_lunar_soi_transit_time(lunar_results, verbose=False)
        elliptical_results = hyperbolic_to_elliptical_conversion(lunar_results, verbose=False)
        
        # Check results
        if elliptical_results.get('error') == 'COLLISION_TRAJECTORY':
            print("❌ COLLISION TRAJECTORY - Parameters not safe!")
            return False
        
        # Calculate totals
        total_delta_v = departure_results['delta_v_departure_kms'] + elliptical_results['total_delta_v']
        total_time = geo_results['tof_hours'] + soi_transit_results['soi_transit_time_hours']
        
        print("✅ SUCCESS - No collision detected!")
        print(f"Total Mission ΔV: {total_delta_v:.3f} km/s")
        print(f"  Earth Departure: {departure_results['delta_v_departure_kms']:.3f} km/s")
        print(f"  Lunar Operations: {elliptical_results['total_delta_v']:.3f} km/s")
        print(f"Total Mission Time: {total_time:.1f} hours")
        print(f"Final Orbit: {elliptical_results.get('orbit_type', 'Unknown')}")
        
        return True
        
    except Exception as e:
        print(f"❌ ERROR: {e}")
        return False

if __name__ == "__main__":
    success = test_safe_parameters()
    if success:
        print("\n✅ Safe parameters confirmed - Option 3 should work correctly!")
    else:
        print("\n❌ Parameters need adjustment")