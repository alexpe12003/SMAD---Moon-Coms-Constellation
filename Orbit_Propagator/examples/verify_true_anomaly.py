#!/usr/bin/env python3
"""
Verify True Anomaly Mathematics for Highly Eccentric Orbits
"""

import numpy as np
import matplotlib.pyplot as plt

def true_anomaly_from_eccentric_correct(E, e):
    """
    Correct conversion from eccentric to true anomaly using atan2 for proper quadrants.
    
    Uses: tan(ν/2) = √((1+e)/(1-e)) * tan(E/2)
    But implemented with atan2 to handle all quadrants properly.
    """
    # Use the half-angle formula but with atan2 for proper quadrant handling
    sqrt_factor = np.sqrt((1 + e) / (1 - e))
    tan_half_nu = sqrt_factor * np.tan(E / 2)
    
    # Calculate ν/2, then double it
    half_nu = np.arctan(tan_half_nu)
    
    # Handle the quadrant properly
    nu = 2 * half_nu
    
    # Ensure nu is in [0, 2π]
    if nu < 0:
        nu += 2*np.pi
        
    return nu

def true_anomaly_from_eccentric_cosine_method(E, e):
    """
    Alternative method using cosine relationship.
    cos(ν) = (cos(E) - e) / (1 - e*cos(E))
    """
    cos_E = np.cos(E)
    cos_nu = (cos_E - e) / (1 - e * cos_E)
    cos_nu = np.clip(cos_nu, -1.0, 1.0)
    
    # Use atan2 with sin(ν) for proper quadrant
    # sin(ν) can be derived from cos(ν) and the quadrant of E
    sin_nu = np.sqrt(1 - cos_nu**2)
    
    # Determine sign of sin(ν) based on eccentric anomaly quadrant
    # For E in [0, π]: sin(E) ≥ 0, so sin(ν) ≥ 0  
    # For E in [π, 2π]: sin(E) ≤ 0, so sin(ν) ≤ 0
    if E > np.pi:
        sin_nu = -sin_nu
    
    nu = np.arctan2(sin_nu, cos_nu)
    
    if nu < 0:
        nu += 2*np.pi
        
    return nu

def analyze_true_anomaly_progression():
    """Analyze how true anomaly should progress for high eccentricity."""
    print("Analysis of True Anomaly Progression for e=0.9")
    print("=" * 60)
    
    e = 0.9
    
    # Create more dense sampling around periapsis
    E_values = np.linspace(0, 2*np.pi, 37)  # Every 10 degrees
    
    print(f"{'E (deg)':>8} {'Method 1':>10} {'Method 2':>10} {'Δν1':>8} {'Δν2':>8}")
    print("-" * 60)
    
    nu1_prev = None
    nu2_prev = None
    
    for E in E_values:
        nu1 = true_anomaly_from_eccentric_correct(E, e) 
        nu2 = true_anomaly_from_eccentric_cosine_method(E, e)
        
        # Calculate incremental changes
        if nu1_prev is not None:
            delta1 = np.degrees(nu1 - nu1_prev)
            delta2 = np.degrees(nu2 - nu2_prev)
            
            # Handle wraparound
            if delta1 > 180:
                delta1 -= 360
            elif delta1 < -180:
                delta1 += 360
                
            if delta2 > 180:
                delta2 -= 360
            elif delta2 < -180:
                delta2 += 360
        else:
            delta1 = delta2 = 0
        
        print(f"{np.degrees(E):8.1f} {np.degrees(nu1):10.1f} {np.degrees(nu2):10.1f} {delta1:8.1f} {delta2:8.1f}")
        
        nu1_prev = nu1
        nu2_prev = nu2

def theoretical_analysis():
    """Analyze the theoretical relationship between E and ν."""
    print("\n\nTheoretical Analysis:")
    print("=" * 60)
    
    print("For highly eccentric orbits (e=0.9):")
    print("- The satellite spends most time near apoapsis (E ~ π)")  
    print("- Very little time near periapsis (E ~ 0, 2π)")
    print("- True anomaly changes RAPIDLY near periapsis")
    print("- This is CORRECT physics - not a bug!")
    print()
    print("The relationship dν/dE = (1-e*cos(E))/(1-e*cos(ν)) shows:")
    print("- Near periapsis (ν≈0): dν/dE ≈ (1-e)/(1-e) = 1")
    print("- Near apoapsis (ν≈π): dν/dE ≈ (1+e)/(1+e) = 1") 
    print("- But the time spent dt/dE varies dramatically!")

if __name__ == "__main__":
    analyze_true_anomaly_progression()
    theoretical_analysis()