#!/usr/bin/env python3
"""
Test script to identify true anomaly calculation issues
"""

import sys
import os
import numpy as np

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.Kepler import true_anomaly_from_eccentric
from core.RK4 import calculate_true_anomaly_from_state
from core.Constants import MU_EARTH

def test_true_anomaly_conversion():
    """Test the true anomaly conversion function."""
    print("Testing true anomaly conversion from eccentric anomaly:")
    print("=" * 60)
    
    e = 0.9  # High eccentricity like Orbit B
    
    # Test several eccentric anomalies around the orbit
    E_values = np.array([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi, 
                        5*np.pi/4, 3*np.pi/2, 7*np.pi/4, 2*np.pi])
    
    print(f"Eccentricity: {e}")
    print(f"{'E (rad)':>8} {'E (deg)':>8} {'ν (rad)':>8} {'ν (deg)':>8} {'Δν (deg)':>10}")
    print("-" * 60)
    
    nu_prev = None
    for E in E_values:
        nu = true_anomaly_from_eccentric(E, e)
        delta_nu = 0 if nu_prev is None else np.degrees(nu - nu_prev)
        print(f"{E:8.3f} {np.degrees(E):8.1f} {nu:8.3f} {np.degrees(nu):8.1f} {delta_nu:10.1f}")
        nu_prev = nu
    
    print("\nLooking for discontinuous jumps in true anomaly...")
    
def test_state_vector_true_anomaly():
    """Test true anomaly calculation from state vectors."""
    print("\nTesting true anomaly from state vectors:")
    print("=" * 60)
    
    # Simple circular orbit test
    r = 7000  # km
    v = np.sqrt(MU_EARTH / r)  # circular velocity
    
    # Test positions around the orbit
    angles = np.linspace(0, 2*np.pi, 9)  # 0 to 360 degrees in 45-degree steps
    
    print(f"{'Angle (deg)':>12} {'ν calc (deg)':>12} {'Error (deg)':>12}")
    print("-" * 40)
    
    for angle in angles:
        # Position and velocity for circular orbit
        pos = r * np.array([np.cos(angle), np.sin(angle), 0])
        vel = v * np.array([-np.sin(angle), np.cos(angle), 0])
        
        nu_calc = calculate_true_anomaly_from_state(pos, vel, MU_EARTH)
        error = np.degrees(nu_calc - angle)
        
        print(f"{np.degrees(angle):12.1f} {np.degrees(nu_calc):12.1f} {error:12.3f}")

def analyze_orbit_b_true_anomaly():
    """Analyze true anomaly behavior for Orbit B parameters."""
    print("\nAnalyzing Orbit B true anomaly behavior:")
    print("=" * 60)
    
    # Orbit B parameters
    a = 70000  # km
    e = 0.9
    
    # Test eccentric anomalies that should give smooth true anomaly progression
    E_test = np.linspace(0, 2*np.pi, 13)  # Every 30 degrees in E
    
    print(f"Semimajor axis: {a} km")
    print(f"Eccentricity: {e}")
    print(f"{'E (deg)':>8} {'ν (deg)':>8} {'Δν (deg)':>10}")
    print("-" * 30)
    
    nu_prev = None
    for E in E_test:
        nu = true_anomaly_from_eccentric(E, e)
        
        # Calculate change in true anomaly
        if nu_prev is not None:
            delta_nu = nu - nu_prev
            # Handle 2π wraparound
            if delta_nu > np.pi:
                delta_nu -= 2*np.pi
            elif delta_nu < -np.pi:
                delta_nu += 2*np.pi
        else:
            delta_nu = 0
            
        print(f"{np.degrees(E):8.1f} {np.degrees(nu):8.1f} {np.degrees(delta_nu):10.1f}")
        nu_prev = nu

if __name__ == "__main__":
    test_true_anomaly_conversion()
    test_state_vector_true_anomaly()
    analyze_orbit_b_true_anomaly()