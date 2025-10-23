#!/usr/bin/env python3
"""
Debug the lunar_results to see what fields are available
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from operations.earth_operations import calculate_earth_departure_delta_v
from operations.trajectory_calculations import lunar_trajectory_calculations
from operations.lunar_operations import lunar_soi_calculations, calculate_lunar_soi_transit_time, hyperbolic_to_elliptical_conversion

def debug_lunar_results():
    """
    Debug the lunar_results structure
    """
    print("Debugging lunar_results structure...")
    print("=" * 50)
    
    # Known working parameters
    R0 = 1.050
    V0 = 1.372
    gamma0 = 0.0
    lambda1 = 240.0
    
    print(f"Testing with: γ0={gamma0}°, λ1={lambda1}°")
    
    try:
        # Calculate up to lunar results
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
        
        print("\nlunar_results keys:")
        for key in sorted(lunar_results.keys()):
            value = lunar_results[key]
            print(f"  {key}: {value}")
        
        # Test different gamma0 values to see collision behavior
        print(f"\nTesting different gamma0 values:")
        
        for test_gamma0 in [0.0, 5.0, 10.0, 15.0]:
            try:
                test_geo_results = lunar_trajectory_calculations(R0, V0, test_gamma0, lambda1, verbose=False)
                test_lunar_results = lunar_soi_calculations(
                    test_geo_results['r1'], 
                    test_geo_results['v1'], 
                    test_geo_results['phi1_deg'],
                    lambda1,
                    test_geo_results['gamma1_deg'],
                    verbose=False
                )
                
                # Check for perigee radius
                rp_hyp = None
                for key in test_lunar_results.keys():
                    if 'rp' in key.lower() or 'perigee' in key.lower():
                        rp_hyp = test_lunar_results[key]
                        print(f"  γ0={test_gamma0}°: {key} = {rp_hyp}")
                        break
                
                if rp_hyp is None:
                    print(f"  γ0={test_gamma0}°: No perigee radius found!")
                    
            except Exception as e:
                print(f"  γ0={test_gamma0}°: ERROR - {e}")
                
    except Exception as e:
        print(f"ERROR: {e}")

if __name__ == "__main__":
    debug_lunar_results()