"""
Lunar-related calculations for transfer missions
"""

import math
from config import *


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


def hyperbolic_to_elliptical_conversion(lunar_results, target_perigee_altitude_km=600, verbose=True):
    """
    Calculate the maneuver required to convert from hyperbolic flyby to elliptical orbit
    
    Parameters:
    - lunar_results: Dictionary containing lunar trajectory analysis results
    - target_perigee_altitude_km: Target perigee altitude in km
    - verbose: Whether to print detailed output
    
    Returns:
    - Dictionary containing conversion analysis results
    """
    
    # Extract lunar trajectory parameters
    rp_hyp = lunar_results['rp']  # Current hyperbolic perigee radius
    e_hyp = lunar_results['e_lunar']  # Current hyperbolic eccentricity
    h_hyp = lunar_results['h_lunar']  # Current angular momentum
    hp_hyp = lunar_results['hp']  # Current hyperbolic perigee altitude
    
    # Target elliptical orbit parameters
    rp_ell = MOON_RADIUS_KM + target_perigee_altitude_km  # Target perigee radius
    ra_ell = rp_hyp  # Target apoapsis = current perigee
    ha_ell = hp_hyp  # Target apoapsis altitude
    
    # Elliptical orbit parameters
    a_ell = (rp_ell + ra_ell) / 2  # Semi-major axis
    e_ell = (ra_ell - rp_ell) / (ra_ell + rp_ell)  # Eccentricity
    p_ell = a_ell * (1 - e_ell**2)  # Semi-latus rectum
    h_ell = math.sqrt(MU_MOON_KMS * p_ell)  # Angular momentum
    
    # Velocities in elliptical orbit
    v_ell_perigee = math.sqrt(MU_MOON_KMS * (2/rp_ell - 1/a_ell))  # At perigee
    v_ell_apoapsis = math.sqrt(MU_MOON_KMS * (2/ra_ell - 1/a_ell))  # At apoapsis
    
    # Velocity in original hyperbolic trajectory at perigee
    v_hyp_perigee = math.sqrt((1 + e_hyp) * MU_MOON_KMS / rp_hyp)
    
    # Calculate required delta-V at apoapsis (current perigee location)
    delta_v_conversion = v_hyp_perigee - v_ell_apoapsis  # CORRECTED: Should be v_hyp - v_ell, not abs()
    
    # Orbital period and transfer time
    T_ell = 2 * math.pi * math.sqrt(a_ell**3 / MU_MOON_KMS)  # Orbital period in seconds
    t_perigee_to_apoapsis = T_ell / 2  # Time from perigee to apoapsis
    
    # Circularization at perigee
    v_circular_perigee = math.sqrt(MU_MOON_KMS / rp_ell)
    delta_v_circularization = v_ell_perigee - v_circular_perigee  # CORRECTED: Should be v_ell - v_circular
    
    # Total delta-V (use absolute values for total)
    total_delta_v = abs(delta_v_conversion) + abs(delta_v_circularization)
    
    if verbose:
        print("=" * 60)
        print("HYPERBOLIC TO ELLIPTICAL TRAJECTORY CONVERSION")
        print("=" * 60)
        print(f"Current hyperbolic trajectory:")
        print(f"  Perigee radius = {rp_hyp:.0f} km")
        print(f"  Perigee altitude = {hp_hyp:.0f} km above surface")
        print(f"  Eccentricity = {e_hyp:.3f}")
        print(f"  Angular momentum = {h_hyp:.0f} km²/sec")
        print(f"  Hyperbolic excess velocity = {lunar_results['v2_kms']:.3f} km/sec")
        
        print(f"\nTarget elliptical trajectory:")
        print(f"  Perigee radius = {rp_ell:.0f} km")
        print(f"  Perigee altitude = {target_perigee_altitude_km:.0f} km above surface")
        print(f"  Apoapsis radius = {ra_ell:.0f} km")
        print(f"  Apoapsis altitude = {ha_ell:.0f} km above surface")
        
        print(f"\nElliptical orbit parameters:")
        print(f"  Semi-major axis = {a_ell:.0f} km")
        print(f"  Eccentricity = {e_ell:.3f}")
        print(f"  Semi-latus rectum = {p_ell:.0f} km")
        print(f"  Angular momentum = {h_ell:.0f} km²/sec")
        
        print(f"\nElliptical orbit velocities:")
        print(f"  Velocity at perigee = {v_ell_perigee:.3f} km/sec")
        print(f"  Velocity at apoapsis = {v_ell_apoapsis:.3f} km/sec")
        
        print(f"\nOriginal hyperbolic trajectory:")
        print(f"  Velocity at perigee = √(1 + e) √(μm/rp) = {v_hyp_perigee:.3f} km/sec")
        
        print(f"\n--- TRAJECTORY CONVERSION MANEUVER ---")
        print(f"Maneuver location: Current perigee (future apoapsis)")
        print(f"Altitude of maneuver: {hp_hyp:.0f} km above Moon surface")
        print(f"Delta-V required: {abs(delta_v_conversion):.3f} km/sec")
        
        if delta_v_conversion < 0:
            print(f"Maneuver type: Retrograde burn (deceleration)")
        else:
            print(f"Maneuver type: Prograde burn (acceleration)")
        
        print(f"\nElliptical orbit characteristics:")
        print(f"  Orbital period = {T_ell:.0f} seconds = {T_ell/3600:.2f} hours")
        print(f"  Time from perigee to apoapsis = {t_perigee_to_apoapsis:.0f} seconds = {t_perigee_to_apoapsis/3600:.2f} hours")
        
        print(f"\n--- CIRCULARIZATION AT PERIGEE ({target_perigee_altitude_km:.0f} km) ---")
        print(f"Circular velocity at {target_perigee_altitude_km:.0f} km altitude = {v_circular_perigee:.3f} km/sec")
        print(f"Elliptical velocity at perigee = {v_ell_perigee:.3f} km/sec")
        print(f"Delta-V for circularization = {abs(delta_v_circularization):.3f} km/sec")
        print(f"Maneuver type: {'Retrograde' if delta_v_circularization > 0 else 'Prograde'} burn ({'deceleration' if delta_v_circularization > 0 else 'acceleration'})")
        
        print(f"\n--- TOTAL MISSION DELTA-V ---")
        print(f"Delta-V for hyperbolic to elliptical = {abs(delta_v_conversion):.3f} km/sec")
        print(f"Delta-V for elliptical to circular = {abs(delta_v_circularization):.3f} km/sec")
        print(f"TOTAL DELTA-V REQUIRED = {total_delta_v:.3f} km/sec")
    
    return {
        'rp_hyperbolic': rp_hyp,
        'hp_hyperbolic': hp_hyp,
        'e_hyperbolic': e_hyp,
        'v_hyp_perigee': v_hyp_perigee,
        'rp_elliptical': rp_ell,
        'ra_elliptical': ra_ell,
        'a_elliptical': a_ell,
        'e_elliptical': e_ell,
        'p_elliptical': p_ell,
        'h_elliptical': h_ell,
        'v_ell_perigee': v_ell_perigee,
        'v_ell_apoapsis': v_ell_apoapsis,
        'delta_v_conversion': delta_v_conversion,
        'delta_v_magnitude': abs(delta_v_conversion),  # Added for compatibility
        'orbital_period_sec': T_ell,
        'orbital_period_hours': T_ell/3600,
        'maneuver_altitude_km': hp_hyp,  # Added for compatibility
        'v_circular_perigee': v_circular_perigee,
        'v_circular_600km': v_circular_perigee,  # Added alias for compatibility
        'delta_v_circularization': delta_v_circularization,
        'delta_v_circ_magnitude': abs(delta_v_circularization),  # Added for compatibility
        'total_delta_v': total_delta_v,
        'target_perigee_altitude_km': target_perigee_altitude_km
    }