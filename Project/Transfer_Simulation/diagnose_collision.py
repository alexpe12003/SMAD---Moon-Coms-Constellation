#!/usr/bin/env python3
"""
Diagnostic script to check the optimal parameters issue
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from operations.earth_operations import calculate_earth_departure_delta_v
from operations.trajectory_calculations import lunar_trajectory_calculations  
from operations.lunar_operations import lunar_soi_calculations

def diagnose_optimal_parameters():
    """Check what's wrong with the optimal parameters"""
    print("DIAGNOSING OPTIMAL PARAMETERS ISSUE")
    print("=" * 50)
    
    # Optimal parameters from the image
    R0 = 1.050      # DU
    V0 = 1.372      # DU/TU
    gamma0 = 29.5   # degrees
    lambda1 = 36.4  # degrees
    
    print(f"Testing parameters:")
    print(f"  R0 = {R0} DU")
    print(f"  V0 = {V0} DU/TU") 
    print(f"  γ0 = {gamma0}°")
    print(f"  λ1 = {lambda1}°")
    print()
    
    try:
        # Step 1: Geocentric trajectory
        print("Step 1: Geocentric trajectory calculations...")
        geo_results = lunar_trajectory_calculations(R0, V0, gamma0, lambda1, verbose=False)
        print(f"✓ r1 = {geo_results['r1']:.3f} DU")
        print(f"✓ v1 = {geo_results['v1']:.4f} DU/TU")
        print(f"✓ phi1 = {geo_results['phi1_deg']:.2f}°")
        print(f"✓ gamma1 = {geo_results['gamma1_deg']:.2f}°")
        
        # Step 2: Lunar SOI calculations
        print(f"\nStep 2: Lunar SOI calculations...")
        lunar_results = lunar_soi_calculations(
            geo_results['r1'], 
            geo_results['v1'], 
            geo_results['phi1_deg'],
            lambda1,
            geo_results['gamma1_deg'],
            verbose=False
        )
        
        print(f"✓ Perigee radius (rp): {lunar_results['rp']:.1f} km")
        print(f"✓ Perigee altitude (hp): {lunar_results['hp']:.1f} km")
        print(f"✓ Moon radius: 1737 km")
        print(f"✓ Eccentricity: {lunar_results['e_lunar']:.3f}")
        
        # Check for collision
        if lunar_results['rp'] < 1737:
            print(f"\n❌ COLLISION DETECTED!")
            print(f"   Perigee radius ({lunar_results['rp']:.1f} km) < Moon radius (1737 km)")
            print(f"   Perigee altitude = {lunar_results['hp']:.1f} km (NEGATIVE - INSIDE MOON!)")
            
            # This is the problem - we can't have an orbit inside the Moon
            print(f"\n🔍 ROOT CAUSE ANALYSIS:")
            print(f"   The spacecraft trajectory passes through the Moon's center")
            print(f"   This means the parameters result in a collision course")
            print(f"   We need different parameters that result in a safe flyby")
            
        else:
            print(f"\n✅ Safe flyby - perigee above Moon surface")
            
        # Step 3: Try orbit conversion anyway to see what happens
        print(f"\nStep 3: Attempting orbit conversion...")
        from operations.lunar_operations import hyperbolic_to_elliptical_conversion
        
        if lunar_results['rp'] > 1737:
            elliptical_results = hyperbolic_to_elliptical_conversion(lunar_results, verbose=False)
            print(f"✓ Conversion successful")
            print(f"  Final orbit type: {elliptical_results.get('orbit_type', 'unknown')}")
            print(f"  Delta-V required: {elliptical_results['total_delta_v']:.3f} km/s")
        else:
            print(f"❌ Cannot convert collision trajectory to orbit")
            print(f"   Physics violation: Cannot orbit inside a solid body")
            
    except Exception as e:
        print(f"❌ Error during calculations: {e}")
        import traceback
        traceback.print_exc()

def suggest_fix():
    """Suggest how to fix this issue"""
    print(f"\n" + "=" * 50)
    print("SUGGESTED FIXES:")
    print("=" * 50)
    print("1. The 'optimal' parameters are actually invalid - they cause a collision")
    print("2. We need to find parameters that result in a safe flyby (perigee > 1737 km)")
    print("3. Options:")
    print("   a) Use a working parameter set (like λ1=229.1° from previous tests)")
    print("   b) Modify the optimization to include collision avoidance constraints")
    print("   c) Check if the optimization results were misread or calculated incorrectly")
    print("\n4. Let's test a known safe parameter set for comparison...")
    
    # Test known good parameters
    print(f"\nTesting known safe parameters (λ1=229.1°):")
    try:
        geo_results = lunar_trajectory_calculations(1.05, 1.372, 10.0, 229.1, verbose=False)
        lunar_results = lunar_soi_calculations(
            geo_results['r1'], geo_results['v1'], geo_results['phi1_deg'],
            229.1, geo_results['gamma1_deg'], verbose=False
        )
        print(f"✓ Safe perigee altitude: {lunar_results['hp']:.0f} km")
    except Exception as e:
        print(f"❌ Error with safe parameters: {e}")

if __name__ == "__main__":
    diagnose_optimal_parameters()
    suggest_fix()