#!/usr/bin/env python3
"""
Test option 3 with the corrected display
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

def test_option3_display():
    """Test the corrected display for option 3"""
    print("TESTING OPTION 3 WITH CORRECTED DISPLAY")
    print("=" * 50)
    
    try:
        # Import the function
        from main import perform_optimal_parameters_calculation
        
        print("Running perform_optimal_parameters_calculation()...")
        print("=" * 50)
        
        # Run the function
        results = perform_optimal_parameters_calculation()
        
        if results:
            user_inputs, departure_results, geo_results, lunar_results, soi_transit_results, circular_results = results
            
            print(f"\n🔍 CHECKING ORBIT CONVERSION RESULTS:")
            print(f"  orbit_type: {circular_results.get('orbit_type', 'NOT FOUND')}")
            print(f"  Natural perigee altitude: {lunar_results['hp']:.1f} km")
            print(f"  Expected: elliptical orbit (since {lunar_results['hp']:.1f} < 1500)")
            
            if circular_results.get('orbit_type') == 'elliptical':
                print(f"  ✅ Correct orbit type!")
                print(f"  Final perigee: {circular_results.get('final_hp', 'NOT FOUND')} km")
                print(f"  Final apoapsis: {circular_results.get('final_ha', 'NOT FOUND')} km")
            else:
                print(f"  ❌ Wrong orbit type: {circular_results.get('orbit_type')}")
                
        else:
            print("❌ Function returned None")
            
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_option3_display()