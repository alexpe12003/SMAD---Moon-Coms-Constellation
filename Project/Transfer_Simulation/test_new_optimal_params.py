#!/usr/bin/env python3
"""
Test the new optimal parameters to ensure they don't cause collision
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from operations.trajectory_calculations import lunar_trajectory_calculations
from operations.lunar_operations import lunar_soi_calculations, hyperbolic_to_elliptical_conversion

def test_new_optimal_parameters():
    """Test the new optimal parameters"""
    print("TESTING NEW OPTIMAL PARAMETERS")
    print("=" * 50)
    
    # New optimal parameters from the image
    R0 = 1.050      # DU
    V0 = 1.372      # DU/TU
    gamma0 = 10.5   # degrees
    lambda1 = 236.4 # degrees
    
    print(f"New optimal parameters:")
    print(f"  R0 = {R0} DU")
    print(f"  V0 = {V0} DU/TU")
    print(f"  γ0 = {gamma0}°")
    print(f"  λ1 = {lambda1}°")
    
    try:
        # Test trajectory calculations
        print(f"\nStep 1: Geocentric trajectory...")
        geo_results = lunar_trajectory_calculations(R0, V0, gamma0, lambda1, verbose=False)
        print(f"✓ Geocentric trajectory successful")
        
        print(f"Step 2: Lunar SOI calculations...")
        lunar_results = lunar_soi_calculations(
            geo_results['r1'], geo_results['v1'], geo_results['phi1_deg'],
            lambda1, geo_results['gamma1_deg'], verbose=False
        )
        
        perigee_radius = lunar_results['rp']
        perigee_altitude = lunar_results['hp']
        
        print(f"✓ Lunar SOI calculations successful")
        print(f"  Perigee radius: {perigee_radius:.1f} km")
        print(f"  Perigee altitude: {perigee_altitude:.1f} km")
        print(f"  Moon radius: 1737 km")
        
        # Check for collision
        if perigee_radius < 1737:
            print(f"❌ COLLISION WARNING: Perigee inside Moon!")
            print(f"   These parameters are unsafe!")
        else:
            print(f"✅ Safe trajectory: Perigee above Moon surface")
        
        print(f"Step 3: Orbit conversion...")
        orbit_results = hyperbolic_to_elliptical_conversion(lunar_results, verbose=False)
        
        if orbit_results.get('error') == 'COLLISION_TRAJECTORY':
            print(f"❌ Orbit conversion detected collision")
        else:
            print(f"✅ Orbit conversion successful")
            print(f"  Orbit type: {orbit_results.get('orbit_type', 'unknown')}")
            print(f"  Total Delta-V: {orbit_results['total_delta_v']:.3f} km/s")
            
            if orbit_results.get('orbit_type') == 'elliptical':
                print(f"  Final orbit: {orbit_results['final_hp']:.0f} km × {orbit_results['final_ha']:.0f} km")
            else:
                print(f"  Final orbit: {orbit_results.get('final_hp', orbit_results['hp_circular']):.0f} km circular")
        
        # Comparison with expected logic
        print(f"\n🔍 LOGIC CHECK:")
        if perigee_altitude < 1500:
            print(f"✓ Expected: Elliptical orbit (perigee {perigee_altitude:.1f} < 1500 km)")
            if orbit_results.get('orbit_type') == 'elliptical' and orbit_results.get('final_ha') == 1500:
                print(f"✅ CORRECT: Got elliptical with apoapsis at 1500 km")
            else:
                print(f"❌ INCORRECT: Expected elliptical with apoapsis at 1500 km")
        else:
            print(f"✓ Expected: Circular orbit (perigee {perigee_altitude:.1f} ≥ 1500 km)")
            if orbit_results.get('orbit_type') == 'circular':
                print(f"✅ CORRECT: Got circular orbit")
            else:
                print(f"❌ INCORRECT: Expected circular orbit")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_new_optimal_parameters()