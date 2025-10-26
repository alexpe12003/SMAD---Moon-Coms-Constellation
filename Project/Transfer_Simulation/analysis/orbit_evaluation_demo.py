import numpy as np
import sys
import os

# Add paths for Lambert solver modules
sys.path.append(r'c:\Users\eduar\OneDrive\Anexos de email\Ambiente de Trabalho\SA\PythonForMongoDB\SMAD---Moon-Coms-Constellation\lambert_solvers\duda_implementation\Duda_Lambert\lambert')
sys.path.append(r'c:\Users\eduar\OneDrive\Anexos de email\Ambiente de Trabalho\SA\PythonForMongoDB\SMAD---Moon-Coms-Constellation\lambert_solvers\duda_implementation\Duda_Lambert\Auxiliary')

try:
    from lambert import lambert
    from coe_sv import coe_from_sv
    print("✅ Successfully imported Lambert solver modules\n")
except ImportError as e:
    print(f"❌ Import error: {e}")
    exit(1)

# Constants
mu = 398600.4418  # km^3/s^2, Earth's gravitational parameter

# True anomalies (degrees) for all cases
theta1_deg = 30.0
theta2_deg = 60.0
theta1 = np.radians(theta1_deg)
theta2 = np.radians(theta2_deg)

print(f"Analysis Parameters:")
print(f"θ₁ = {theta1_deg}° = {theta1:.4f} rad")
print(f"θ₂ = {theta2_deg}° = {theta2:.4f} rad")
print(f"μ = {mu:,.1f} km³/s²")
print("="*80)

# Define orbit parameters
orbits = {
    'A': {'e': 0.1, 'a': 7000.0, 'type': 'Elliptical'},     # e < 1, a > 0
    'B': {'e': 0.9, 'a': 70000.0, 'type': 'Elliptical'},   # e < 1, a > 0  
    'C': {'e': 2.0, 'a': -7000.0, 'type': 'Hyperbolic'}    # e > 1, a < 0
}

def calculate_orbit_properties_detailed(e, a, theta1, theta2, orbit_name):
    """Calculate orbital properties with detailed step-by-step output"""
    
    print(f"\n🔧 DETAILED CALCULATIONS FOR ORBIT {orbit_name}")
    print("-" * 60)
    print(f"Input parameters:")
    print(f"  • Eccentricity (e) = {e}")
    print(f"  • Semi-major axis (a) = {a:,.0f} km")
    print(f"  • True anomaly 1 (θ₁) = {np.degrees(theta1):.1f}°")
    print(f"  • True anomaly 2 (θ₂) = {np.degrees(theta2):.1f}°")
    
    print(f"\nStep 1: Calculate radial distances using orbital equation")
    print(f"r = a(1-e²)/(1+e*cos(θ))")
    
    # Calculate r1
    numerator1 = a * (1 - e**2)
    denominator1 = 1 + e * np.cos(theta1)
    r1 = numerator1 / denominator1
    print(f"r₁ = {a}*(1-{e}²)/(1+{e}*cos({np.degrees(theta1):.1f}°))")
    print(f"r₁ = {a}*{1-e**2:.6f}/{denominator1:.6f} = {r1:.2f} km")
    
    # Calculate r2
    numerator2 = a * (1 - e**2)
    denominator2 = 1 + e * np.cos(theta2)
    r2 = numerator2 / denominator2
    print(f"r₂ = {a}*(1-{e}²)/(1+{e}*cos({np.degrees(theta2):.1f}°))")
    print(f"r₂ = {a}*{1-e**2:.6f}/{denominator2:.6f} = {r2:.2f} km")
    
    print(f"\nStep 2: Calculate velocities using vis-viva equation")
    print(f"v = √(μ(2/r - 1/a))")
    
    # Calculate v1
    term1_v1 = 2/r1
    term2_v1 = 1/a
    v1_squared = mu * (term1_v1 - term2_v1)
    v1 = np.sqrt(v1_squared)
    print(f"v₁ = √({mu:.1f}*(2/{r1:.2f} - 1/{a:.0f}))")
    print(f"v₁ = √({mu:.1f}*({term1_v1:.8f} - {term2_v1:.8f}))")
    print(f"v₁ = √{v1_squared:.2f} = {v1:.4f} km/s")
    
    # Calculate v2
    term1_v2 = 2/r2
    term2_v2 = 1/a
    v2_squared = mu * (term1_v2 - term2_v2)
    v2 = np.sqrt(v2_squared)
    print(f"v₂ = √({mu:.1f}*(2/{r2:.2f} - 1/{a:.0f}))")
    print(f"v₂ = √({mu:.1f}*({term1_v2:.8f} - {term2_v2:.8f}))")
    print(f"v₂ = √{v2_squared:.2f} = {v2:.4f} km/s")
    
    print(f"\nStep 3: Calculate time of flight")
    
    # Time of flight calculation
    if e < 1:  # Elliptical orbit
        print(f"Elliptical orbit (e < 1): Using Kepler's equation")
        
        # Eccentric anomaly calculation
        print(f"Calculate eccentric anomalies using: cos(E) = (e + cos(θ))/(1 + e*cos(θ))")
        
        cos_E1 = (e + np.cos(theta1)) / (1 + e * np.cos(theta1))
        E1 = np.arccos(np.clip(cos_E1, -1, 1))
        if theta1 > np.pi:
            E1 = 2 * np.pi - E1
        print(f"cos(E₁) = ({e} + cos({np.degrees(theta1):.1f}°))/(1 + {e}*cos({np.degrees(theta1):.1f}°)) = {cos_E1:.6f}")
        print(f"E₁ = {np.degrees(E1):.4f}° = {E1:.6f} rad")

        cos_E2 = (e + np.cos(theta2)) / (1 + e * np.cos(theta2))
        E2 = np.arccos(np.clip(cos_E2, -1, 1))
        if theta2 > np.pi:
            E2 = 2 * np.pi - E2
        print(f"cos(E₂) = ({e} + cos({np.degrees(theta2):.1f}°))/(1 + {e}*cos({np.degrees(theta2):.1f}°)) = {cos_E2:.6f}")
        print(f"E₂ = {np.degrees(E2):.4f}° = {E2:.6f} rad")

        # Mean anomalies
        print(f"\nCalculate mean anomalies using: M = E - e*sin(E)")
        M1 = E1 - e * np.sin(E1)
        M2 = E2 - e * np.sin(E2)
        print(f"M₁ = {E1:.6f} - {e}*sin({E1:.6f}) = {M1:.6f} rad")
        print(f"M₂ = {E2:.6f} - {e}*sin({E2:.6f}) = {M2:.6f} rad")

        # Period and time of flight
        print(f"\nCalculate orbital period: T = 2π√(a³/μ)")
        T = 2 * np.pi * np.sqrt(a**3 / mu)
        print(f"T = 2π√({a}³/{mu:.1f}) = {T:.2f} s = {T/3600:.4f} hours")
        
        print(f"\nCalculate time of flight: Δt = (M₂ - M₁) * T / (2π)")
        tof = (M2 - M1) * T / (2 * np.pi)
        if tof < 0:
            tof += T
            print(f"Δt was negative, added one period: Δt = {tof:.2f} s")
        print(f"Δt = ({M2:.6f} - {M1:.6f}) * {T:.2f} / (2π) = {tof:.2f} s = {tof/3600:.4f} hours")
            
    else:  # Hyperbolic orbit (e > 1)
        print(f"Hyperbolic orbit (e > 1): Using hyperbolic Kepler's equation")
        
        # Hyperbolic eccentric anomaly
        print(f"Calculate hyperbolic eccentric anomalies using: cosh(F) = (e + cos(θ))/(1 + e*cos(θ))")
        
        cosh_F1 = (e + np.cos(theta1)) / (1 + e * np.cos(theta1))
        F1 = np.arccosh(cosh_F1)
        print(f"cosh(F₁) = ({e} + cos({np.degrees(theta1):.1f}°))/(1 + {e}*cos({np.degrees(theta1):.1f}°)) = {cosh_F1:.6f}")
        print(f"F₁ = {F1:.6f} rad")
        
        cosh_F2 = (e + np.cos(theta2)) / (1 + e * np.cos(theta2))
        F2 = np.arccosh(cosh_F2)
        print(f"cosh(F₂) = ({e} + cos({np.degrees(theta2):.1f}°))/(1 + {e}*cos({np.degrees(theta2):.1f}°)) = {cosh_F2:.6f}")
        print(f"F₂ = {F2:.6f} rad")
        
        # Hyperbolic mean anomalies
        print(f"\nCalculate hyperbolic mean anomalies using: M = e*sinh(F) - F")
        M1 = e * np.sinh(F1) - F1
        M2 = e * np.sinh(F2) - F2
        print(f"M₁ = {e}*sinh({F1:.6f}) - {F1:.6f} = {M1:.6f}")
        print(f"M₂ = {e}*sinh({F2:.6f}) - {F2:.6f} = {M2:.6f}")
        
        # Time of flight for hyperbolic orbit
        print(f"\nCalculate mean motion: n = √(μ/|a|³)")
        n = np.sqrt(mu / (-a)**3)  # Mean motion (note: a is negative)
        print(f"n = √({mu:.1f}/{-a}³) = {n:.8f} rad/s")
        
        print(f"\nCalculate time of flight: Δt = (M₂ - M₁) / n")
        tof = (M2 - M1) / n
        print(f"Δt = ({M2:.6f} - {M1:.6f}) / {n:.8f} = {tof:.2f} s = {tof/3600:.4f} hours")
    
    print(f"\n📊 FINAL RESULTS:")
    print(f"  • r₁ = {r1:.2f} km")
    print(f"  • r₂ = {r2:.2f} km")
    print(f"  • v₁ = {v1:.4f} km/s")
    print(f"  • v₂ = {v2:.4f} km/s")
    print(f"  • Time of flight = {tof:.2f} s = {tof/3600:.4f} hours")
    
    return r1, r2, v1, v2, tof

def solve_lambert_problem_detailed(r1, r2, theta1, theta2, tof, orbit_name):
    """Solve Lambert's problem with detailed output"""
    
    print(f"\n🎯 LAMBERT SOLVER FOR ORBIT {orbit_name}")
    print("-" * 60)
    
    # Construct position vectors
    print(f"Step 1: Construct position vectors")
    R1 = np.array([r1 * np.cos(theta1), r1 * np.sin(theta1), 0])
    R2 = np.array([r2 * np.cos(theta2), r2 * np.sin(theta2), 0])
    print(f"R₁ = [{R1[0]:.2f}, {R1[1]:.2f}, {R1[2]:.2f}] km")
    print(f"R₂ = [{R2[0]:.2f}, {R2[1]:.2f}, {R2[2]:.2f}] km")
    
    print(f"\nStep 2: Solve Lambert's problem")
    print(f"Input to Lambert solver:")
    print(f"  • R₁ = {R1}")
    print(f"  • R₂ = {R2}")
    print(f"  • Time of flight = {tof:.2f} s")
    print(f"  • μ = {mu} km³/s²")
    
    try:
        # Solve Lambert problem
        V1_lambert, V2_lambert, theta1_lambert, theta2_lambert = lambert(R1, R2, tof, string='pro', mu=mu)
        
        print(f"\nStep 3: Lambert solution results")
        print(f"V₁_Lambert = [{V1_lambert[0]:.4f}, {V1_lambert[1]:.4f}, {V1_lambert[2]:.4f}] km/s")
        print(f"V₂_Lambert = [{V2_lambert[0]:.4f}, {V2_lambert[1]:.4f}, {V2_lambert[2]:.4f}] km/s")
        print(f"|V₁_Lambert| = {np.linalg.norm(V1_lambert):.4f} km/s")
        print(f"|V₂_Lambert| = {np.linalg.norm(V2_lambert):.4f} km/s")
        print(f"θ₁_Lambert = {np.degrees(theta1_lambert):.4f}°")
        print(f"θ₂_Lambert = {np.degrees(theta2_lambert):.4f}°")
        
        # Calculate orbital elements from Lambert solution
        print(f"\nStep 4: Calculate transfer orbit elements")
        coe_lambert = coe_from_sv(R1, V1_lambert, mu)
        print(f"Transfer orbit eccentricity = {coe_lambert[1]:.6f}")
        print(f"Transfer orbit semi-major axis = {coe_lambert[6]:.2f} km")
        
        return {
            'success': True,
            'V1_magnitude': np.linalg.norm(V1_lambert),
            'V2_magnitude': np.linalg.norm(V2_lambert),
            'theta1_lambert': np.degrees(theta1_lambert),
            'theta2_lambert': np.degrees(theta2_lambert),
            'lambert_e': coe_lambert[1],
            'lambert_a': coe_lambert[6],
            'V1': V1_lambert,
            'V2': V2_lambert
        }
    except Exception as error:
        print(f"❌ Lambert solver failed: {error}")
        return {'success': False, 'error': str(error)}

# Main Analysis
print("DETAILED ORBIT EVALUATION PROCESS")
print("="*80)

for orbit_name, params in orbits.items():
    print(f"\n\n{'='*20} ORBIT {orbit_name} - {params['type'].upper()} {'='*20}")
    
    e = params['e']
    a = params['a']
    
    # Calculate basic orbital properties with detailed steps
    r1, r2, v1, v2, tof = calculate_orbit_properties_detailed(e, a, theta1, theta2, orbit_name)
    
    # Solve Lambert problem with detailed steps
    lambert_result = solve_lambert_problem_detailed(r1, r2, theta1, theta2, tof, orbit_name)
    
    if lambert_result['success']:
        print(f"\n🔍 COMPARISON:")
        v1_diff = abs(v1 - lambert_result['V1_magnitude'])
        v2_diff = abs(v2 - lambert_result['V2_magnitude'])
        print(f"Original v₁ = {v1:.4f} km/s vs Lambert v₁ = {lambert_result['V1_magnitude']:.4f} km/s")
        print(f"Difference |v₁ - v₁_Lambert| = {v1_diff:.6f} km/s")
        print(f"Original v₂ = {v2:.4f} km/s vs Lambert v₂ = {lambert_result['V2_magnitude']:.4f} km/s")
        print(f"Difference |v₂ - v₂_Lambert| = {v2_diff:.6f} km/s")

print(f"\n\n{'='*80}")
print("EVALUATION COMPLETE")
print("="*80)