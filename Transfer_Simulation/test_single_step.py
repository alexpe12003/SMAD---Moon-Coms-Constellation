#!/usr/bin/env python3
"""
Simple test of the optimization components with known good values
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from operations.earth_operations import calculate_earth_departure_delta_v
from operations.trajectory_calculations import lunar_trajectory_calculations
from operations.lunar_operations import lunar_soi_calculations, calculate_lunar_soi_transit_time, hyperbolic_to_elliptical_conversion

def test_single_optimization_step():
    """Test a single optimization step with known good values"""
    print("Testing single optimization step with known good values...")
    print("=" * 60)
    
    # Known good values from our previous tests
    R0 = 1.05
    V0 = 1.372
    gamma0 = 10.0
    lambda1 = 229.1
    
    try:
        print(f"Testing: R0={R0}, V0={V0}, γ0={gamma0}°, λ1={lambda1}°")
        
        # Earth departure analysis
        departure_results = calculate_earth_departure_delta_v(R0, V0, verbose=False)
        earth_departure_dv = departure_results['delta_v_departure_kms']
        print(f"✓ Earth departure ΔV: {earth_departure_dv:.3f} km/s")
        
        # Geocentric trajectory calculations
        geo_results = lunar_trajectory_calculations(R0, V0, gamma0, lambda1, verbose=False)
        print(f"✓ Geocentric trajectory TOF: {geo_results['tof_hours']:.1f} hours")
        
        # Lunar SOI calculations
        lunar_results = lunar_soi_calculations(
            geo_results['r1'], 
            geo_results['v1'], 
            geo_results['phi1_deg'],
            lambda1,
            geo_results['gamma1_deg'],
            verbose=False
        )
        print(f"✓ Lunar SOI perigee altitude: {lunar_results['hp']:.0f} km")
        
        # SOI transit time
        soi_transit_results = calculate_lunar_soi_transit_time(lunar_results, verbose=False)
        print(f"✓ SOI transit time: {soi_transit_results['soi_transit_time_hours']:.2f} hours")
        
        # Orbit conversion (new function)
        elliptical_results = hyperbolic_to_elliptical_conversion(lunar_results, verbose=False)
        lunar_dv = elliptical_results['total_delta_v']
        print(f"✓ Lunar orbit insertion ΔV: {lunar_dv:.3f} km/s")
        print(f"✓ Orbit type: {elliptical_results.get('orbit_type', 'unknown')}")
        
        # Total mission delta-V
        total_delta_v = earth_departure_dv + lunar_dv
        total_time = geo_results['tof_hours'] + soi_transit_results['soi_transit_time_hours']
        
        print(f"\n🎯 RESULTS:")
        print(f"   Earth departure ΔV: {earth_departure_dv:.3f} km/s")
        print(f"   Lunar insertion ΔV: {lunar_dv:.3f} km/s") 
        print(f"   TOTAL MISSION ΔV: {total_delta_v:.3f} km/s")
        print(f"   Total mission time: {total_time:.1f} hours")
        
        return True
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_single_optimization_step()
    if success:
        print("\n✅ Single step test successful!")
    else:
        print("\n❌ Single step test failed!")