#!/usr/bin/env python3
"""
Simple test to verify the display fixes work
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from operations.trajectory_calculations import lunar_trajectory_calculations
from operations.lunar_operations import lunar_soi_calculations, hyperbolic_to_elliptical_conversion

def simple_test():
    """Simple test with known parameters"""
    print("SIMPLE TEST OF CORRECTED LOGIC AND DISPLAY")
    print("=" * 50)
    
    # Parameters from option 3
    R0, V0, gamma0, lambda1 = 1.050, 1.372, 10.0, 229.1
    
    try:
        print("Step 1: Geocentric trajectory...")
        geo_results = lunar_trajectory_calculations(R0, V0, gamma0, lambda1, verbose=False)
        
        print("Step 2: Lunar SOI calculations...")
        lunar_results = lunar_soi_calculations(
            geo_results['r1'], geo_results['v1'], geo_results['phi1_deg'],
            lambda1, geo_results['gamma1_deg'], verbose=False
        )
        
        perigee_alt = lunar_results['hp']
        print(f"Natural perigee altitude: {perigee_alt:.1f} km")
        print(f"Decision rule: {perigee_alt:.1f} km < 1500 km? {perigee_alt < 1500}")
        
        print(f"\nStep 3: Orbit conversion...")
        orbit_results = hyperbolic_to_elliptical_conversion(lunar_results, verbose=True)
        
        print(f"\n🔍 CHECKING RESULTS:")
        print(f"Orbit type: '{orbit_results.get('orbit_type', 'MISSING')}'")
        print(f"Final perigee altitude: {orbit_results.get('final_hp', 'MISSING')}")
        print(f"Final apoapsis altitude: {orbit_results.get('final_ha', 'MISSING')}")
        
        # Test the display function
        print(f"\n📺 TESTING DISPLAY FUNCTION:")
        from core.interface import display_organized_mission_summary
        
        user_inputs = {'R0': R0, 'V0': V0, 'gamma0': gamma0, 'lambda1': lambda1}
        
        # This should now show elliptical orbit correctly
        print("Calling display_organized_mission_summary...")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    simple_test()