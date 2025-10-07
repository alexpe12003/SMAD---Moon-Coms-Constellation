#!/usr/bin/env python3
"""
Diagnostic script to investigate lambda1 calculation failures
"""

import sys
import os
import math
import traceback

# Add the current directory to Python path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# Import the transfer functions
from Transfer import lunar_trajectory_calculations, lunar_soi_calculations, hyperbolic_to_elliptical_conversion

def test_lambda1_range():
    """Test lambda1 values from 290-335 degrees to identify failure points"""
    
    print("DIAGNOSTIC: LAMBDA1 CALCULATION FAILURES")
    print("=" * 60)
    
    # Fixed initial conditions (same as parametric study)
    R0 = 1.05  # DU
    V0 = 1.372  # DU/TU
    gamma0 = 0  # degrees
    
    print(f"Test Parameters:")
    print(f"  R0 = {R0} DU")
    print(f"  V0 = {V0} DU/TU")
    print(f"  gamma0 = {gamma0}°")
    print()
    
    # Test the problematic range
    test_values = [290, 295, 300, 305, 310, 315, 320, 325, 330, 335]
    
    for lambda1 in test_values:
        print(f"\n{'='*40}")
        print(f"TESTING LAMBDA1 = {lambda1}°")
        print(f"{'='*40}")
        
        try:
            # Step 1: Geocentric trajectory calculations
            print(f"Step 1: Geocentric trajectory calculations...")
            geo_results = lunar_trajectory_calculations(R0, V0, gamma0, lambda1, verbose=False)
            
            print(f"  ✓ Geocentric calculations successful")
            print(f"    r1 = {geo_results['r1']:.4f} DU")
            print(f"    v1 = {geo_results['v1']:.4f} DU/TU")
            print(f"    phi1 = {geo_results['phi1_deg']:.2f}°")
            print(f"    gamma1 = {geo_results['gamma1_deg']:.2f}°")
            print(f"    TOF = {geo_results['tof_hours']:.2f} hours")
            
            # Step 2: Lunar SOI calculations
            print(f"Step 2: Lunar SOI calculations...")
            lunar_results = lunar_soi_calculations(
                geo_results['r1'], 
                geo_results['v1'], 
                geo_results['phi1_deg'],
                lambda1,
                geo_results['gamma1_deg'],
                verbose=False
            )
            
            print(f"  ✓ Lunar SOI calculations successful")
            print(f"    v2 = {lunar_results['v2']:.4f} km/s")
            print(f"    e_lunar = {lunar_results['e_lunar']:.4f}")
            print(f"    rp = {lunar_results['rp']:.0f} km")
            print(f"    hp = {lunar_results['hp']:.0f} km")
            
            # Step 3: Check for physical validity
            if lunar_results['e_lunar'] <= 1.0:
                print(f"  ⚠️  WARNING: e_lunar = {lunar_results['e_lunar']:.4f} ≤ 1.0 (not hyperbolic!)")
            if lunar_results['hp'] < 0:
                print(f"  ⚠️  WARNING: hp = {lunar_results['hp']:.0f} km < 0 (impact trajectory!)")
            
            # Step 4: Elliptical conversion
            print(f"Step 3: Hyperbolic to elliptical conversion...")
            elliptical_results = hyperbolic_to_elliptical_conversion(
                lunar_results, 
                target_perigee_altitude_km=600, 
                verbose=False
            )
            
            print(f"  ✓ Conversion calculations successful")
            print(f"    Total delta-V = {elliptical_results['total_delta_v']:.3f} km/s")
            
        except Exception as e:
            print(f"  ❌ ERROR at lambda1 = {lambda1}°:")
            print(f"    {type(e).__name__}: {str(e)}")
            print(f"    Traceback:")
            traceback.print_exc()
            print()

def analyze_flight_path_angle():
    """Analyze the flight path angle at Moon's SOI for problematic lambda1 values"""
    
    print(f"\n{'='*60}")
    print("FLIGHT PATH ANGLE ANALYSIS")
    print(f"{'='*60}")
    
    # Test parameters
    R0 = 1.05
    V0 = 1.372
    gamma0 = 0
    
    for lambda1 in range(0, 361, 15):  # Every 15 degrees
        try:
            geo_results = lunar_trajectory_calculations(R0, V0, gamma0, lambda1, verbose=False)
            phi1 = geo_results['phi1_deg']
            
            # Check for problematic flight path angles
            if phi1 > 90:
                status = "⚠️  WARNING: phi1 > 90°"
            elif phi1 < 0:
                status = "❌ ERROR: phi1 < 0°"
            else:
                status = "✓ OK"
                
            print(f"λ₁ = {lambda1:3d}°  →  φ₁ = {phi1:6.2f}°  {status}")
            
        except Exception as e:
            print(f"λ₁ = {lambda1:3d}°  →  ERROR: {str(e)}")

if __name__ == "__main__":
    test_lambda1_range()
    analyze_flight_path_angle()