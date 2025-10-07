#!/usr/bin/env python3
"""
Detailed analysis of the math domain error in lunar_soi_calculations
"""

import math

def analyze_asin_domain_error():
    """Analyze the asin domain error for problematic lambda1 values"""
    
    print("DETAILED ANALYSIS: math.asin() DOMAIN ERROR")
    print("=" * 60)
    
    # Test parameters from the failing cases
    test_cases = [
        {"lambda1": 290, "v1_kms": 0.1140 * 7.905, "gamma1_deg": -9.77, "phi1_deg": 77.32},
        {"lambda1": 295, "v1_kms": 0.1165 * 7.905, "gamma1_deg": -9.57, "phi1_deg": 77.39},
        {"lambda1": 300, "v1_kms": 0.1189 * 7.905, "gamma1_deg": -9.28, "phi1_deg": 77.46},
        {"lambda1": 305, "v1_kms": 0.1212 * 7.905, "gamma1_deg": -8.91, "phi1_deg": 77.51},
        {"lambda1": 310, "v1_kms": 0.1235 * 7.905, "gamma1_deg": -8.45, "phi1_deg": 77.56},
        {"lambda1": 315, "v1_kms": 0.1256 * 7.905, "gamma1_deg": -7.91, "phi1_deg": 77.61},
        {"lambda1": 320, "v1_kms": 0.1277 * 7.905, "gamma1_deg": -7.28, "phi1_deg": 77.64},
        {"lambda1": 325, "v1_kms": 0.1295 * 7.905, "gamma1_deg": -6.57, "phi1_deg": 77.67},
        {"lambda1": 330, "v1_kms": 0.1312 * 7.905, "gamma1_deg": -5.79, "phi1_deg": 77.70},
    ]
    
    # Constants
    vm = 1.018  # Moon's orbital velocity in km/sec
    
    print("The problematic equation is:")
    print("ε2 = arcsin[(vm/v2)*cos(λ1) - (v1/v2)*cos(λ1 + γ1 - φ1)]")
    print()
    print("For arcsin() to work, the argument must be in range [-1, 1]")
    print()
    
    for case in test_cases:
        lambda1 = case["lambda1"]
        v1_kms = case["v1_kms"]
        gamma1_deg = case["gamma1_deg"]
        phi1_deg = case["phi1_deg"]
        
        # Convert to radians
        lambda1_rad = math.radians(lambda1)
        gamma1_rad = math.radians(gamma1_deg)
        phi1_rad = math.radians(phi1_deg)
        
        # Calculate v2 (relative velocity to Moon)
        v2 = math.sqrt(v1_kms**2 + vm**2 - 2*v1_kms*vm*math.cos(phi1_rad))
        
        # Calculate the argument for arcsin
        term1 = (vm/v2) * math.cos(lambda1_rad)
        term2 = (v1_kms/v2) * math.cos(lambda1_rad + gamma1_rad - phi1_rad)
        asin_arg = term1 - term2
        
        # Check if it's in valid range
        is_valid = -1.0 <= asin_arg <= 1.0
        status = "✓ Valid" if is_valid else "❌ INVALID"
        
        print(f"λ₁ = {lambda1:3d}°:")
        print(f"  v1 = {v1_kms:.4f} km/s")
        print(f"  v2 = {v2:.4f} km/s")
        print(f"  vm/v2 = {vm/v2:.6f}")
        print(f"  v1/v2 = {v1_kms/v2:.6f}")
        print(f"  cos(λ₁) = {math.cos(lambda1_rad):8.6f}")
        print(f"  cos(λ₁+γ₁-φ₁) = {math.cos(lambda1_rad + gamma1_rad - phi1_rad):8.6f}")
        print(f"  term1 = (vm/v2)*cos(λ₁) = {term1:8.6f}")
        print(f"  term2 = (v1/v2)*cos(λ₁+γ₁-φ₁) = {term2:8.6f}")
        print(f"  asin_arg = term1 - term2 = {asin_arg:8.6f}  {status}")
        if not is_valid:
            print(f"  ERROR: asin({asin_arg:.6f}) is outside [-1, 1] range!")
        print()

def suggest_solution():
    """Suggest solutions for the domain error"""
    
    print("SUGGESTED SOLUTIONS:")
    print("=" * 30)
    print("1. CLAMP METHOD: Limit asin argument to [-1, 1] range")
    print("   asin_arg = max(-1.0, min(1.0, asin_arg))")
    print()
    print("2. SKIP METHOD: Skip invalid lambda1 values in parametric study")
    print("   if not (-1.0 <= asin_arg <= 1.0): continue")
    print()
    print("3. PHYSICS CHECK: These may represent physically impossible")
    print("   geometries where the spacecraft cannot reach the Moon's SOI")
    print("   with the given transfer parameters.")
    print()
    print("The error suggests that for lambda1 values ~295-330°, the")
    print("calculated flight path geometry results in impossible relative")
    print("velocity vectors at the Moon's SOI.")

if __name__ == "__main__":
    analyze_asin_domain_error()
    suggest_solution()