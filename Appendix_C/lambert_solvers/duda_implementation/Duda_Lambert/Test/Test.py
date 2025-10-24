import numpy as np

# Import Lambert solver and auxiliary functions
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../lambert')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../Auxiliary')))
from lambert import lambert
from coe_sv import coe_from_sv

# Constants
mu = 398600.4418  # km^3/s^2, Earth's gravitational parameter

# True anomalies (degrees) for all cases
theta1_deg = 30.0
theta2_deg = 60.0
theta1 = np.radians(theta1_deg)
theta2 = np.radians(theta2_deg)

# Define orbit parameters
orbits = {
    'A': {'e': 0.1, 'a': 7000.0, 'type': 'Elliptical'},     # e < 1, a > 0
    'B': {'e': 0.9, 'a': 70000.0, 'type': 'Elliptical'},   # e < 1, a > 0  
    'C': {'e': 2.0, 'a': -7000.0, 'type': 'Hyperbolic'}    # e > 1, a < 0
}

def calculate_orbit_properties(e, a, theta1, theta2):
    """Calculate orbital properties for given parameters"""
    
    # Basic orbital calculations
    r1 = a * (1 - e**2) / (1 + e * np.cos(theta1))
    r2 = a * (1 - e**2) / (1 + e * np.cos(theta2))
    
    # Velocities using vis-viva equation
    v1 = np.sqrt(mu * (2/r1 - 1/a))
    v2 = np.sqrt(mu * (2/r2 - 1/a))
    
    # Time of flight calculation
    if e < 1:  # Elliptical orbit
        # Eccentric anomaly calculation
        cos_E1 = (e + np.cos(theta1)) / (1 + e * np.cos(theta1))
        E1 = np.arccos(np.clip(cos_E1, -1, 1))
        if theta1 > np.pi:
            E1 = 2 * np.pi - E1

        cos_E2 = (e + np.cos(theta2)) / (1 + e * np.cos(theta2))
        E2 = np.arccos(np.clip(cos_E2, -1, 1))
        if theta2 > np.pi:
            E2 = 2 * np.pi - E2

        # Mean anomalies
        M1 = E1 - e * np.sin(E1)
        M2 = E2 - e * np.sin(E2)

        # Period and time of flight
        T = 2 * np.pi * np.sqrt(a**3 / mu)
        tof = (M2 - M1) * T / (2 * np.pi)
        if tof < 0:
            tof += T
            
    else:  # Hyperbolic orbit (e > 1)
        # Hyperbolic eccentric anomaly
        cosh_F1 = (e + np.cos(theta1)) / (1 + e * np.cos(theta1))
        F1 = np.arccosh(cosh_F1)
        
        cosh_F2 = (e + np.cos(theta2)) / (1 + e * np.cos(theta2))
        F2 = np.arccosh(cosh_F2)
        
        # Hyperbolic mean anomalies
        M1 = e * np.sinh(F1) - F1
        M2 = e * np.sinh(F2) - F2
        
        # Time of flight for hyperbolic orbit
        n = np.sqrt(mu / (-a)**3)  # Mean motion (note: a is negative)
        tof = (M2 - M1) / n
    
    return r1, r2, v1, v2, tof

def solve_lambert_problem(r1, r2, theta1, theta2, tof):
    """Solve Lambert's problem and return results"""
    
    # Construct position vectors
    R1 = np.array([r1 * np.cos(theta1), r1 * np.sin(theta1), 0])
    R2 = np.array([r2 * np.cos(theta2), r2 * np.sin(theta2), 0])
    
    try:
        # Solve Lambert problem
        V1_lambert, V2_lambert, theta1_lambert, theta2_lambert = lambert(R1, R2, tof, string='pro', mu=mu)
        
        # Calculate orbital elements from Lambert solution
        coe_lambert = coe_from_sv(R1, V1_lambert, mu)
        
        return {
            'success': True,
            'V1_magnitude': np.linalg.norm(V1_lambert),
            'V2_magnitude': np.linalg.norm(V2_lambert),
            'theta1_lambert': np.degrees(theta1_lambert),
            'theta2_lambert': np.degrees(theta2_lambert),
            'lambert_e': coe_lambert[1],
            'lambert_a': coe_lambert[6],
            'R1': R1,
            'R2': R2,
            'V1': V1_lambert,
            'V2': V2_lambert
        }
    except Exception as error:
        return {'success': False, 'error': str(error)}

# Main Analysis
print("="*80)
print("ORBITAL ANALYSIS: THREE ORBIT TYPES")
print("="*80)
print(f"Analysis Case: θ₁ = {theta1_deg}°, θ₂ = {theta2_deg}°")
print("="*80)

results = {}

for orbit_name, params in orbits.items():
    print(f"\n🚀 ORBIT {orbit_name} - {params['type'].upper()} ORBIT")
    print("-" * 60)
    
    e = params['e']
    a = params['a']
    
    print(f"📋 Orbital Parameters:")
    print(f"   • Eccentricity (e): {e}")
    print(f"   • Semi-major axis (a): {a:,.0f} km")
    print(f"   • Orbit type: {params['type']}")
    
    try:
        # Calculate basic orbital properties
        r1, r2, v1, v2, tof = calculate_orbit_properties(e, a, theta1, theta2)
        
        print(f"\n📊 Calculated Properties:")
        print(f"   • Initial radius (r₁): {r1:,.2f} km")
        print(f"   • Final radius (r₂): {r2:,.2f} km") 
        print(f"   • Initial velocity (v₁): {v1:.2f} km/s")
        print(f"   • Final velocity (v₂): {v2:.2f} km/s")
        print(f"   • Time of flight: {tof:,.2f} s ({tof/3600:.2f} hours)")
        
        # Solve Lambert problem
        lambert_result = solve_lambert_problem(r1, r2, theta1, theta2, tof)
        
        if lambert_result['success']:
            print(f"\n🎯 Lambert Solution:")
            print(f"   • Initial velocity magnitude: {lambert_result['V1_magnitude']:.2f} km/s")
            print(f"   • Final velocity magnitude: {lambert_result['V2_magnitude']:.2f} km/s")
            print(f"   • Lambert θ₁: {lambert_result['theta1_lambert']:.2f}°")
            print(f"   • Lambert θ₂: {lambert_result['theta2_lambert']:.2f}°")
            print(f"   • Transfer orbit eccentricity: {lambert_result['lambert_e']:.4f}")
            print(f"   • Transfer orbit semi-major axis: {lambert_result['lambert_a']:,.2f} km")
            
            # Velocity comparison
            v1_diff = abs(v1 - lambert_result['V1_magnitude'])
            v2_diff = abs(v2 - lambert_result['V2_magnitude'])
            print(f"\n🔍 Velocity Comparison:")
            print(f"   • |v₁ - v₁_Lambert|: {v1_diff:.4f} km/s")
            print(f"   • |v₂ - v₂_Lambert|: {v2_diff:.4f} km/s")
            
        else:
            print(f"\n❌ Lambert Solution Failed: {lambert_result['error']}")
            
        # Store results
        results[orbit_name] = {
            'e': e, 'a': a, 'r1': r1, 'r2': r2, 'v1': v1, 'v2': v2, 'tof': tof,
            'lambert': lambert_result
        }
        
    except Exception as error:
        print(f"\n❌ Error calculating orbit {orbit_name}: {error}")
        results[orbit_name] = {'error': str(error)}

# Summary Table
print("\n" + "="*80)
print("📋 SUMMARY TABLE")
print("="*80)

header = f"{'Orbit':<6} {'Type':<12} {'e':<8} {'a (km)':<12} {'r₁ (km)':<10} {'r₂ (km)':<10} {'ToF (h)':<10} {'Status':<10}"
print(header)
print("-" * 80)

for orbit_name, data in results.items():
    if 'error' not in data:
        orbit_info = orbits[orbit_name]
        status = "✅ OK" if data['lambert']['success'] else "❌ FAIL"
        row = f"{orbit_name:<6} {orbit_info['type']:<12} {data['e']:<8} {data['a']:>10,.0f} {data['r1']:>8,.0f} {data['r2']:>8,.0f} {data['tof']/3600:>8.2f} {status:<10}"
        print(row)
    else:
        print(f"{orbit_name:<6} {'ERROR':<12} {'-':<8} {'-':<12} {'-':<10} {'-':<10} {'-':<10} {'❌ ERR':<10}")

print("\n" + "="*80)