"""
Trajectory calculations for Earth-Moon transfer
"""

import math
from config import *


def lunar_trajectory_calculations(R0, V0, gamma0_deg, lambda1_deg, verbose=True):
    """
    Calculate geocentric trajectory parameters for lunar transfer
    
    Parameters:
    - R0: Initial radius in DU
    - V0: Initial velocity in DU/TU
    - gamma0_deg: Initial flight path angle in degrees
    - lambda1_deg: Lambda1 angle in degrees
    - verbose: Whether to print detailed output
    
    Returns:
    - Dictionary containing trajectory analysis results
    """
    
    # Convert angles to radians
    gamma0_rad = math.radians(gamma0_deg)
    lambda1_rad = math.radians(lambda1_deg)
    
    # Constants
    D = MOON_DISTANCE  # DU
    Rs = MOON_SOI_RADIUS  # DU
    
    if verbose:
        print("=" * 60)
        print("DETAILED CALCULATION STEPS")
        print("=" * 60)
        print(f"Input parameters:")
        print(f"R0 = {R0} DU")
        print(f"V0 = {V0} DU/TU")
        print(f"φ0 (gamma0) = {gamma0_deg}°")
        print(f"λ1 = {lambda1_deg}°")
        print(f"\nConstants:")
        print(f"D = {D} DU")
        print(f"Rs = {Rs} DU")
    
    # Step 1: Calculate specific energy
    epsilon = V0**2 / 2 - MU_EARTH / R0
    
    if verbose:
        print(f"\nStep 1 - Specific Energy (ε):")
        print(f"ε = V0²/2 - 1/R0 = {epsilon:.4f} DU²/TU²")
    
    # Step 2: Calculate specific angular momentum
    h = R0 * V0 * math.cos(gamma0_rad)
    
    if verbose:
        print(f"\nStep 2 - Specific Angular Momentum (h):")
        print(f"h = R0 * V0 * cos(γ0) = {h:.4f} DU²/TU")
    
    # Step 3: Calculate distance at Moon's SOI
    r1 = math.sqrt(D**2 + Rs**2 - 2*D*Rs*math.cos(lambda1_rad))
    
    if verbose:
        print(f"\nStep 3 - Distance at Moon's SOI (r1):")
        print(f"r1 = √(D² + Rs² - 2*D*Rs*cos(λ1)) = {r1:.4f} DU")
    
    # Step 4: Calculate velocity at Moon's SOI
    v1 = math.sqrt(2*(epsilon + h**2/(2*r1**2) + MU_EARTH/r1))
    
    if verbose:
        print(f"\nStep 4 - Velocity at Moon's SOI (v1):")
        print(f"v1 = √(2(ε + h²/(2r1²) + 1/r1)) = {v1:.4f} DU/TU")
    
    # Step 5: Calculate flight path angle at Moon's SOI
    phi1_rad = math.acos(h / (r1 * v1))
    phi1_deg = math.degrees(phi1_rad)
    
    if verbose:
        print(f"\nStep 5 - Flight Path Angle at Moon's SOI (φ1):")
        print(f"φ1 = arccos(h/(r1*v1)) = {phi1_deg:.2f}°")
    
    # Step 6: Calculate phase angle at Moon's SOI
    gamma1_rad = math.asin((Rs/r1) * math.sin(lambda1_rad))
    gamma1_deg = math.degrees(gamma1_rad)
    
    if verbose:
        print(f"\nStep 6 - Phase Angle at Moon's SOI (γ1):")
        print(f"γ1 = arcsin((Rs/r1) * sin(λ1)) = {gamma1_deg:.2f}° = {gamma1_rad:.3f} rad")
    
    # Step 7: Calculate orbital parameters
    p = h**2  # Semi-latus rectum
    a = -1 / (2 * epsilon)  # Semi-major axis
    e = math.sqrt(1 + 2*epsilon*h**2)  # Eccentricity
    
    if verbose:
        print(f"\nStep 7 - Geocentric Trajectory Parameters:")
        print(f"p = h² = {p:.3f} DU")
        print(f"a = -1/(2ε) = {a:.2f} DU")
        print(f"e = √(1 + 2εh²) = {e:.3f}")
    
    # Step 8: Calculate true anomalies
    # Using orbit equation: r = p/(1 + e*cos(v))
    # Solving for true anomaly: cos(v) = (p/r - 1)/e
    
    cos_v0 = (p/R0 - 1) / e
    cos_v1 = (p/r1 - 1) / e
    
    v0_rad = math.acos(cos_v0) if gamma0_deg == 0 else math.acos(cos_v0)
    v1_rad = math.acos(cos_v1)
    
    v0_deg = math.degrees(v0_rad)
    v1_deg = math.degrees(v1_rad)
    
    if verbose:
        print(f"\nStep 8 - Anomaly Calculations:")
        print(f"v0 = {v0_deg:.0f}° (since φ0 = {gamma0_deg}°)")
        print(f"v1 = {v1_deg:.2f}° = {v1_rad:.3f} rad")
    
    # Step 9: Calculate eccentric anomalies
    E0_rad = 2 * math.atan(math.sqrt((1-e)/(1+e)) * math.tan(v0_rad/2))
    E1_rad = 2 * math.atan(math.sqrt((1-e)/(1+e)) * math.tan(v1_rad/2))
    
    E0_deg = math.degrees(E0_rad)
    E1_deg = math.degrees(E1_rad)
    
    if verbose:
        print(f"E0 = {E0_deg:.0f}° (since v0 = {v0_deg:.0f}°)")
        print(f"E1 = {E1_deg:.1f}° = {E1_rad:.2f} rad")
    
    # Step 10: Calculate time of flight using Kepler's equation
    sin_E1 = math.sin(E1_rad)
    sin_E0 = math.sin(E0_rad)
    
    tof_TU = math.sqrt(a**3) * ((E1_rad - e*sin_E1) - (E0_rad - e*sin_E0))
    tof_hours = tof_TU * TU_TO_HOURS
    
    if verbose:
        print(f"\nStep 9 - Time of Flight Calculation:")
        print(f"sin E1 = {sin_E1:.3f}")
        print(f"TOF = √(a³) * [(E1 - e*sin E1) - (E0 - e*sin E0)]")
        print(f"TOF = {tof_TU:.2f} TU = {tof_hours:.2f} hours")
    
    # Step 11: Calculate initial phase angle
    # Using relationship: γ0 = v1 - v0 - γ1 - ωm*(t1 - t0)
    # where ωm is Moon's angular velocity (approximately 2π/month)
    omega_m = 2 * math.pi / (27.32 * 24)  # rad/hour (sidereal month)
    gamma0_calc_rad = v1_rad - v0_rad - gamma1_rad - omega_m * tof_hours
    gamma0_calc_deg = math.degrees(gamma0_calc_rad)
    
    if verbose:
        print(f"\nStep 10 - Initial Phase Angle (γ0):")
        print(f"γ0 = v1 - v0 - γ1 - ωm*(t1 - t0)")
        print(f"γ0 = {gamma0_calc_deg:.2f}° = {gamma0_calc_rad:.3f} rad")
    
    return {
        'R0': R0,
        'V0': V0,
        'gamma0_deg': gamma0_deg,
        'lambda1_deg': lambda1_deg,
        'epsilon': epsilon,
        'h': h,
        'r1': r1,
        'v1': v1,
        'phi1_deg': phi1_deg,
        'phi1_rad': phi1_rad,
        'gamma1_deg': gamma1_deg,
        'gamma1_rad': gamma1_rad,
        'p': p,
        'a': a,
        'e': e,
        'v0_deg': v0_deg,
        'v1_deg': v1_deg,
        'v0_rad': v0_rad,
        'v1_rad': v1_rad,
        'E0_deg': E0_deg,
        'E1_deg': E1_deg,
        'E0_rad': E0_rad,
        'E1_rad': E1_rad,
        'tof_TU': tof_TU,
        'tof_hours': tof_hours,
        'gamma0_calc_deg': gamma0_calc_deg,
        'gamma0_calc_rad': gamma0_calc_rad
    }