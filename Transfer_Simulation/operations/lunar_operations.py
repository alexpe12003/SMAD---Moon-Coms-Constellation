"""
Lunar-related calculations for transfer missions
"""

import math
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import *


def calculate_lunar_soi_transit_time(lunar_results, verbose=True):
    """
    Calculate the time for spacecraft to transit from Moon's SOI entry to perigee
    using hyperbolic orbital mechanics.
    
    Parameters:
    - lunar_results: Dictionary containing lunar trajectory analysis results
    - verbose: Whether to print detailed output
    
    Returns:
    - Dictionary containing SOI transit time analysis
    """
    
    # Extract lunar trajectory parameters
    R_s_km = MOON_SOI_RADIUS * DU_TO_KM  # SOI radius in km
    mu_m_kms = MU_MOON_KMS  # Moon's gravitational parameter in km³/s²
    e_lunar = lunar_results['e_lunar']
    rp = lunar_results['rp']  # periapsis radius in km
    h_lunar = lunar_results['h_lunar']  # angular momentum in km²/sec
    
    # Calculate semi-major axis of hyperbolic trajectory
    a_hyp = rp / (e_lunar - 1)  # For hyperbolic orbits, a is negative, but we use |a|
    
    # True anomaly at SOI entry
    # At SOI boundary: r = R_s = a(e*cosh(F) - 1) for hyperbolic orbit
    # Solving for true anomaly at SOI entry
    cos_nu_soi = (h_lunar**2 / (mu_m_kms * R_s_km) - 1) / e_lunar
    nu_soi_rad = math.acos(cos_nu_soi)
    nu_soi_deg = math.degrees(nu_soi_rad)
    
    # Hyperbolic eccentric anomaly at SOI entry
    F_soi = 2 * math.atanh(math.sqrt((e_lunar - 1)/(e_lunar + 1)) * math.tan(nu_soi_rad/2))
    
    # Hyperbolic eccentric anomaly at perigee (F = 0)
    F_perigee = 0
    
    # Mean motion for hyperbolic orbit
    n_hyp = math.sqrt(mu_m_kms / a_hyp**3)
    
    # Time from SOI entry to perigee using hyperbolic Kepler's equation
    # M = e*sinh(F) - F, and M = n*(t - t0)
    M_soi = e_lunar * math.sinh(F_soi) - F_soi
    M_perigee = e_lunar * math.sinh(F_perigee) - F_perigee  # = 0 at perigee
    
    # Time difference
    delta_M = M_perigee - M_soi  # This will be negative (time before perigee)
    soi_transit_time_sec = abs(delta_M) / n_hyp
    soi_transit_time_hours = soi_transit_time_sec / 3600
    soi_transit_time_minutes = soi_transit_time_sec / 60
    
    if verbose:
        print(f"Hyperbolic Trajectory in Moon's SOI:")
        print(f"  SOI radius: {R_s_km:.0f} km")
        print(f"  Periapsis radius: {rp:.0f} km")
        print(f"  Eccentricity: {e_lunar:.3f}")
        print(f"  Semi-major axis: {a_hyp:.0f} km")
        
        print(f"\nSOI Entry Point:")
        print(f"  True anomaly at SOI: {nu_soi_deg:.1f}°")
        print(f"  Hyperbolic eccentric anomaly: {F_soi:.3f}")
        
        print(f"\nTransit Time Analysis:")
        print(f"  Time from SOI entry to perigee: {soi_transit_time_hours:.2f} hours")
        print(f"  Time from SOI entry to perigee: {soi_transit_time_minutes:.1f} minutes")
        print(f"  Time from SOI entry to perigee: {soi_transit_time_sec:.0f} seconds")
    
    return {
        'soi_radius_km': R_s_km,
        'periapsis_radius_km': rp,
        'eccentricity': e_lunar,
        'semi_major_axis_km': a_hyp,
        'nu_soi_deg': nu_soi_deg,
        'nu_soi_rad': nu_soi_rad,
        'F_soi': F_soi,
        'soi_transit_time_sec': soi_transit_time_sec,
        'soi_transit_time_minutes': soi_transit_time_minutes,
        'soi_transit_time_hours': soi_transit_time_hours
    }


def lunar_soi_calculations(r1, v1, phi1_deg, lambda1_deg, gamma1_deg, verbose=True):
    """
    Calculate lunar sphere of influence trajectory parameters
    
    Parameters:
    - r1: Distance at Moon's SOI in DU
    - v1: Velocity at Moon's SOI in DU/TU
    - phi1_deg: Flight path angle at Moon's SOI in degrees
    - lambda1_deg: Lambda1 angle in degrees
    - gamma1_deg: Gamma1 angle in degrees
    - verbose: Whether to print detailed output
    
    Returns:
    - Dictionary containing lunar SOI analysis results
    """
    
    # Convert units
    v1_kms = v1 * DU_TU_TO_KM_S  # Convert to km/sec
    Rs_km = MOON_SOI_RADIUS * DU_TO_KM  # SOI radius in km
    vm_kms = MOON_ORBITAL_VELOCITY  # Moon's orbital velocity in km/sec
    
    if verbose:
        print("=" * 60)
        print("LUNAR SPHERE OF INFLUENCE CALCULATIONS")
        print("=" * 60)
        print(f"v1 = {v1:.4f} DU/TU = {v1_kms:.3f} km/sec")
        print(f"Rs = {MOON_SOI_RADIUS:.3f} DU = {Rs_km:.0f} km")
        print(f"μm = {MU_MOON_KMS} (Earth canonical units)")
        print(f"vm = {vm_kms:.3f} km/sec")
    
    # Convert angles to radians
    phi1_rad = math.radians(phi1_deg)
    lambda1_rad = math.radians(lambda1_deg)
    gamma1_rad = math.radians(gamma1_deg)
    
    # Calculate relative velocity at Moon's SOI using law of cosines
    # v2² = v1² + vm² - 2*v1*vm*cos(φ1)
    v2_kms = math.sqrt(v1_kms**2 + vm_kms**2 - 2*v1_kms*vm_kms*math.cos(phi1_rad))
    
    # Calculate flight path angle relative to Moon using law of sines
    # ε2 = arcsin[(vm/v2)*cos(λ1) - (v1/v2)*cos(λ1 + γ1 - φ1)]
    angle_term = lambda1_rad + gamma1_rad - phi1_rad
    epsilon2_rad = math.asin((vm_kms/v2_kms)*math.cos(lambda1_rad) - (v1_kms/v2_kms)*math.cos(angle_term))
    epsilon2_deg = math.degrees(epsilon2_rad)
    
    if verbose:
        print(f"\nRelative velocity at Moon's SOI:")
        print(f"v2 = √(v1² + vm² - 2*v1*vm*cos(φ1)) = {v2_kms:.3f} km/sec")
        print(f"ε2 = arcsin[(vm/v2)*cos(λ1) - (v1/v2)*cos(λ1 + γ1 - φ1)] = {epsilon2_deg:.2f}°")
    
    # Lunar-centric trajectory parameters
    h_lunar = Rs_km * v2_kms * math.sin(epsilon2_rad)  # Angular momentum in km²/sec (CORRECTED: should be sin, not cos)
    epsilon_lunar = v2_kms**2 / 2 - MU_MOON_KMS / Rs_km  # Specific energy
    p_lunar = h_lunar**2 / MU_MOON_KMS  # Semi-latus rectum in km
    e_lunar = math.sqrt(1 + 2*epsilon_lunar*h_lunar**2/MU_MOON_KMS**2)  # Eccentricity
    rp_lunar = p_lunar / (1 + e_lunar)  # Periapsis radius in km
    hp_lunar = rp_lunar - MOON_RADIUS_KM  # Periapsis altitude above Moon surface
    
    if verbose:
        print(f"\nLunar-centric trajectory parameters:")
        print(f"h = {h_lunar:.0f} km²/sec")
        print(f"ε = v2²/2 - μm/Rs = {epsilon_lunar:.4f}")
        print(f"p = h²/μm = {p_lunar:.0f} km")
        print(f"e = √(1 + 2εh²/μm²) = {e_lunar:.3f}")
        print(f"rp = p/(1 + e) = {rp_lunar:.0f} km")
        print(f"hp = rp - Rmoon = {hp_lunar:.0f} km")
        print(f"\nThis is the minimum distance of approach to the Moon's surface.")
    
    return {
        'v1_kms': v1_kms,
        'v2': v2_kms,  # Added missing v2 field for compatibility
        'v2_kms': v2_kms,
        'epsilon2_deg': epsilon2_deg,
        'epsilon2_rad': epsilon2_rad,
        'h_lunar': h_lunar,
        'epsilon_lunar': epsilon_lunar,
        'p_lunar': p_lunar,
        'e_lunar': e_lunar,
        'rp': rp_lunar,
        'hp': hp_lunar
    }


def hyperbolic_to_elliptical_conversion(lunar_results, target_perigee_altitude_km=None, verbose=True):
    """
    Convert hyperbolic trajectory to either elliptical or circular orbit
    
    Logic:
    - If natural perigee altitude < 1500km: Create ellipse with natural perigee and apoapsis at 1500km  
    - If natural perigee altitude ≥ 1500km: Circularize at natural perigee
    - Collision detection: Returns error if perigee is inside Moon (< 1737 km radius)
    
    Parameters:
    - lunar_results: Dictionary containing lunar trajectory analysis results
    - target_perigee_altitude_km: Ignored parameter (kept for compatibility)
    - verbose: Whether to print detailed output
    
    Returns:
    - Dictionary containing conversion analysis results, or error dict for collisions
    """
    
    # Extract lunar trajectory parameters
    rp_hyp = lunar_results['rp']  # Current hyperbolic perigee radius
    e_hyp = lunar_results['e_lunar']  # Current hyperbolic eccentricity
    h_hyp = lunar_results['h_lunar']  # Current angular momentum
    hp_hyp = lunar_results['hp']  # Current hyperbolic perigee altitude
    
    # COLLISION DETECTION: Check if trajectory passes through Moon
    if rp_hyp < MOON_RADIUS_KM:
        if verbose:
            print("=" * 60)
            print("❌ COLLISION TRAJECTORY DETECTED!")
            print("=" * 60)
            print(f"Hyperbolic perigee radius: {rp_hyp:.1f} km")
            print(f"Moon radius: {MOON_RADIUS_KM:.1f} km")
            print(f"Perigee altitude: {hp_hyp:.1f} km (NEGATIVE - INSIDE MOON!)")
            print(f"❌ ERROR: Cannot create orbit from collision trajectory")
            print(f"   The spacecraft would crash into the Moon's surface.")
            print(f"   Different mission parameters are required for a safe flyby.")
        
        # Return error result
        return {
            'error': 'COLLISION_TRAJECTORY',
            'error_message': f'Hyperbolic perigee ({rp_hyp:.1f} km) is below Moon surface ({MOON_RADIUS_KM} km)',
            'rp_hyperbolic': rp_hyp,
            'hp_hyperbolic': hp_hyp,
            'collision': True,
            'total_delta_v': float('inf'),  # Infinite delta-V for impossible trajectory
            'orbit_type': 'collision'
        }
    
    # Velocity in original hyperbolic trajectory at perigee
    v_hyp_perigee = math.sqrt((1 + e_hyp) * MU_MOON_KMS / rp_hyp)
    
    # Define target apoapsis altitude
    TARGET_APOAPSIS_ALTITUDE_KM = 1500  # km
    
    if hp_hyp < TARGET_APOAPSIS_ALTITUDE_KM:
        # Mode 1: Create elliptical orbit with apoapsis at 1500km
        # Use natural perigee as perigee, 1500km as apoapsis
        rp_elliptical = rp_hyp  # Current perigee becomes new perigee  
        ra_elliptical = MOON_RADIUS_KM + TARGET_APOAPSIS_ALTITUDE_KM  # Apoapsis at 1500km
        
        # Calculate elliptical orbit parameters
        a_elliptical = (rp_elliptical + ra_elliptical) / 2  # Semi-major axis
        e_elliptical = (ra_elliptical - rp_elliptical) / (ra_elliptical + rp_elliptical)  # Eccentricity
        
        # Calculate velocity at perigee of elliptical orbit
        v_peri_elliptical = math.sqrt(MU_MOON_KMS * (2/rp_elliptical - 1/a_elliptical))
        
        # Delta-V calculation at natural perigee
        delta_v_circularization = v_peri_elliptical - v_hyp_perigee
        total_delta_v = abs(delta_v_circularization)
        
        # Final orbit parameters (elliptical)
        final_rp = rp_elliptical
        final_hp = hp_hyp  # Natural perigee altitude
        final_ra = ra_elliptical
        final_ha = TARGET_APOAPSIS_ALTITUDE_KM
        final_velocity = v_peri_elliptical
        maneuver_altitude = hp_hyp
        orbit_type = "elliptical"
        
    else:
        # Mode 2: Direct circularization at natural perigee (perigee > 1500km)
        v_circular_perigee = math.sqrt(MU_MOON_KMS / rp_hyp)
        delta_v_circularization = v_hyp_perigee - v_circular_perigee  # Deceleration needed
        total_delta_v = abs(delta_v_circularization)
        
        # Final orbit parameters (circular at natural perigee)
        final_rp = rp_hyp
        final_hp = hp_hyp
        final_ra = rp_hyp  # Circular orbit
        final_ha = hp_hyp
        final_velocity = v_circular_perigee
        maneuver_altitude = hp_hyp
        orbit_type = "circular"
    
    if verbose:
        print("=" * 60)
        print("HYPERBOLIC TO ORBIT CONVERSION")
        print("=" * 60)
        print(f"Current hyperbolic trajectory:")
        print(f"  Perigee radius = {rp_hyp:.0f} km")
        print(f"  Perigee altitude = {hp_hyp:.0f} km above surface")
        print(f"  Eccentricity = {e_hyp:.3f}")
        print(f"  Angular momentum = {h_hyp:.0f} km²/sec")
        print(f"  Hyperbolic excess velocity = {lunar_results['v2_kms']:.3f} km/sec")
        
        print(f"\nDecision logic:")
        print(f"  Target apoapsis altitude = {TARGET_APOAPSIS_ALTITUDE_KM} km")
        if hp_hyp < TARGET_APOAPSIS_ALTITUDE_KM:
            print(f"  Natural perigee ({hp_hyp:.0f} km) < {TARGET_APOAPSIS_ALTITUDE_KM} km → Creating elliptical orbit with apoapsis at {TARGET_APOAPSIS_ALTITUDE_KM} km")
        else:
            print(f"  Natural perigee ({hp_hyp:.0f} km) ≥ {TARGET_APOAPSIS_ALTITUDE_KM} km → Circularizing at natural perigee")
        
        print(f"\nVelocities at perigee:")
        print(f"  Hyperbolic velocity = √(1 + e) √(μm/rp) = {v_hyp_perigee:.3f} km/sec")
        print(f"  Target orbit velocity = {final_velocity:.3f} km/sec")
        
        print(f"\n--- ORBITAL MANEUVER ---")
        print(f"Maneuver altitude: {maneuver_altitude:.0f} km above Moon surface")
        print(f"Delta-V required: {total_delta_v:.3f} km/sec")
        if delta_v_circularization < 0:
            print(f"Maneuver type: Retrograde burn (deceleration)")
        else:
            print(f"Maneuver type: Prograde burn (acceleration)")
        
        print(f"\n--- FINAL ORBIT CHARACTERISTICS ({orbit_type.upper()}) ---")
        if orbit_type == "elliptical":
            final_period = 2 * math.pi * math.sqrt(((final_rp + final_ra) / 2)**3 / MU_MOON_KMS)
            print(f"  Perigee radius = {final_rp:.0f} km")
            print(f"  Perigee altitude = {final_hp:.0f} km above surface")
            print(f"  Apoapsis radius = {final_ra:.0f} km")
            print(f"  Apoapsis altitude = {final_ha:.0f} km above surface")
            print(f"  Eccentricity = {(final_ra - final_rp)/(final_ra + final_rp):.3f}")
        else:
            final_period = 2 * math.pi * math.sqrt(final_rp**3 / MU_MOON_KMS)
            print(f"  Orbital radius = {final_rp:.0f} km")
            print(f"  Orbital altitude = {final_hp:.0f} km above surface")
        
        print(f"  Velocity at perigee = {final_velocity:.3f} km/sec")
        print(f"  Orbital period = {final_period:.0f} seconds = {final_period/3600:.2f} hours")
        
        print(f"\n--- TOTAL MISSION DELTA-V ---")
        print(f"TOTAL DELTA-V REQUIRED = {total_delta_v:.3f} km/sec")
    
    # Calculate final orbit period for completeness
    if orbit_type == "elliptical":
        final_period_sec = 2 * math.pi * math.sqrt(((final_rp + final_ra) / 2)**3 / MU_MOON_KMS)
    else:
        final_period_sec = 2 * math.pi * math.sqrt(final_rp**3 / MU_MOON_KMS)
    
    return {
        'rp_hyperbolic': rp_hyp,
        'hp_hyperbolic': hp_hyp,
        'e_hyperbolic': e_hyp,
        'v_hyp_perigee': v_hyp_perigee,
        'final_rp': final_rp,
        'final_hp': final_hp,
        'final_ra': final_ra,
        'final_ha': final_ha,
        'final_velocity': final_velocity,
        'orbit_type': orbit_type,
        'delta_v_circularization': delta_v_circularization,
        'delta_v_magnitude': abs(delta_v_circularization),
        'total_delta_v': total_delta_v,
        'final_period_sec': final_period_sec,
        'final_period_hours': final_period_sec/3600,
        'maneuver_altitude_km': maneuver_altitude,
        'target_apoapsis_altitude_km': TARGET_APOAPSIS_ALTITUDE_KM,
        # Legacy compatibility fields
        'rp_circular': final_rp,
        'hp_circular': final_hp,
        'v_circular': final_velocity,
        'circular_period_sec': final_period_sec,
        'circular_period_hours': final_period_sec/3600,
        'circularization_altitude_km': final_hp
    }