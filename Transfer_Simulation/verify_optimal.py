#!/usr/bin/env python3
"""
Simple verification of the new optimal parameters function
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

def test_optimal_parameters():
    """Test the optimal parameters calculation directly"""
    print("Testing optimal parameters calculation...")
    print("=" * 50)
    
    try:
        # Import the function
        from main import perform_optimal_parameters_calculation
        
        print("✓ Function imported successfully")
        
        # Test the function (but catch any issues)
        print("Running optimal parameters calculation...")
        results = perform_optimal_parameters_calculation()
        
        if results:
            print("✅ Optimal parameters calculation completed successfully!")
            user_inputs, departure_results, geo_results, lunar_results, soi_transit_results, circular_results = results
            
            # Show key results
            total_dv = departure_results['delta_v_departure_kms'] + circular_results['total_delta_v']
            print(f"\n🎯 KEY RESULTS:")
            print(f"   Total Mission Delta-V: {total_dv:.3f} km/s")
            print(f"   Earth-Moon Transit: {geo_results['tof_hours']:.1f} hours")
            print(f"   Natural Perigee: {lunar_results['hp']:.0f} km")
            print(f"   Final Orbit: {circular_results.get('orbit_type', 'unknown')}")
            
        else:
            print("❌ Function returned None")
            
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_optimal_parameters()