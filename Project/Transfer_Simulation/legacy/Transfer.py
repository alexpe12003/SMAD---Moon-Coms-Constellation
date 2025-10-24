import math
import matplotlib.pyplot as plt
import numpy as np

def get_user_input():
    """
    Get mission parameters including lambda1
    Returns a dictionary with the input values
    """
    print("Please enter the following parameters:")
    print("-" * 50)
    
    variables = {}
    
    # Set predefined mission parameters
    variables['R0'] = 1.05  # DU - Both parking orbit radius AND transfer trajectory perigee
    variables['V0'] = 1.372  # DU/TU - Transfer trajectory velocity at perigee
    variables['gamma0'] = 0  # degrees - Flight path angle at perigee
    
    print(f"Using predefined mission parameters:")
    print(f"Parking orbit radius = Transfer perigee (R0) = {variables['R0']} DU")
    print(f"Parking orbit altitude = {(variables['R0']-1)*6378:.0f} km")
    print(f"Transfer trajectory velocity (V0) = {variables['V0']} DU/TU")
    print(f"Transfer flight path angle (gamma0) = {variables['gamma0']}°")
    print()
    
    # Only ask for lambda1 from user
    while True:
        try:
            value = float(input("Enter lambda1 (Point where geocentric trajectory crosses lunar sphere of influence in degrees): "))
            variables['lambda1'] = value
            break
        except ValueError:
            print("Error: Please enter a valid number.")
            continue
    
    # Commented out original user input code
    # # Define the variables to collect
    # var_names = ['R0', 'V0', 'gamma0', 'lambda1']
    # var_descriptions = [
    #     'R0 (Transfer trajectory perigee in DU)',
    #     'V0 (Transfer trajectory velocity in DU/TU)', 
    #     'gamma0 (Transfer flight path angle in degrees)',
    #     'lambda1 (Point where geocentric trajectory crosses lunar sphere of influence in degrees)'
    # ]
    # 
    # # Collect each variable with input validation
    # for var_name, description in zip(var_names, var_descriptions):
    #     while True:
    #         try:
    #             value = float(input(f"Enter {description}: "))
    #             variables[var_name] = value
    #             break
    #         except ValueError:
    #             print("Error: Please enter a valid number.")
    #             continue
    
    return variables

def calculate_earth_departure_delta_v(R0, V0, verbose=True):
    """
    Calculate the delta-V required to depart from Earth circular parking orbit
    to the transfer trajectory. The parking orbit radius equals the transfer perigee (R0).
    
    Parameters:
    R0: Parking orbit radius AND transfer trajectory perigee radius in DU
    V0: Transfer trajectory velocity at perigee in DU/TU
    verbose: whether to print detailed calculations
    
    Returns:
    Dictionary with departure parameters and delta-V
    """
    if verbose:
        print("\n" + "=" * 60)
        print("EARTH DEPARTURE ANALYSIS")
        print("=" * 60)
    
    # Calculate circular velocity at parking orbit (same as transfer perigee)
    v_circular = math.sqrt(1 / R0)  # In canonical units, μ_earth = 1
    
    # The transfer velocity at this radius is given as V0
    v_transfer_at_departure = V0
    
    # Delta-V for departure
    delta_v_departure = v_transfer_at_departure - v_circular
    
    # Convert to km/s (1 DU/TU ≈ 7.905 km/s)
    v_circular_kms = v_circular * 7.905
    v_transfer_kms = v_transfer_at_departure * 7.905
    delta_v_departure_kms = delta_v_departure * 7.905
    
    # Calculate transfer trajectory parameters for reference
    epsilon_transfer = (V0**2 / 2) - (1 / R0)
    a_transfer = -1 / (2 * epsilon_transfer) if epsilon_transfer < 0 else float('inf')
    
    if verbose:
        print(f"Parking Orbit Analysis:")
        print(f"  Radius: {R0:.4f} DU ({(R0-1)*6378:.0f} km altitude)")
        print(f"  Circular velocity: {v_circular:.4f} DU/TU = {v_circular_kms:.3f} km/s")
        
        print(f"\nTransfer Trajectory Analysis:")
        print(f"  Perigee radius: {R0:.4f} DU (same as parking orbit)")
        if epsilon_transfer < 0:
            print(f"  Semi-major axis: {a_transfer:.2f} DU")
            print(f"  Trajectory type: Elliptical")
        else:
            print(f"  Trajectory type: Hyperbolic")
        print(f"  Specific energy: {epsilon_transfer:.4f} DU²/TU²")
        print(f"  Velocity at departure: {v_transfer_at_departure:.4f} DU/TU = {v_transfer_kms:.3f} km/s")
        
        print(f"\nDeparture Maneuver:")
        print(f"  Location: {(R0-1)*6378:.0f} km altitude")
        print(f"  Delta-V required: {delta_v_departure:.4f} DU/TU = {delta_v_departure_kms:.3f} km/s")
        print(f"  Maneuver type: Prograde burn (acceleration)")
    
    return {
        'R0': R0,
        'v_circular': v_circular,
        'v_circular_kms': v_circular_kms,
        'a_transfer': a_transfer if epsilon_transfer < 0 else None,
        'epsilon_transfer': epsilon_transfer,
        'v_transfer_at_departure': v_transfer_at_departure,
        'v_transfer_kms': v_transfer_kms,
        'delta_v_departure': delta_v_departure,
        'delta_v_departure_kms': delta_v_departure_kms,
        'altitude_km': (R0-1)*6378
    }

def calculate_lunar_soi_transit_time(lunar_results, verbose=True):
    """
    Calculate time spent in Moon's SOI from entry to perigee passage
    
    Parameters:
    lunar_results: dictionary with lunar trajectory parameters
    verbose: whether to print detailed calculations
    
    Returns:
    Dictionary with SOI transit time information
    """
    if verbose:
        print("\n" + "=" * 60)
        print("LUNAR SOI TRANSIT TIME ANALYSIS")
        print("=" * 60)
    
    # Moon parameters
    R_s_km = 66300  # Moon's SOI radius in km
    mu_m_kms = 4.903e3  # Moon's gravitational parameter in km³/sec²
    
    # Extract hyperbolic trajectory parameters
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
        'a_hyperbolic': a_hyp,
        'nu_soi_deg': nu_soi_deg,
        'F_soi': F_soi,
        'soi_transit_time_sec': soi_transit_time_sec,
        'soi_transit_time_minutes': soi_transit_time_minutes,
        'soi_transit_time_hours': soi_transit_time_hours
    }

def lunar_trajectory_calculations(R0, V0, gamma0_deg, lambda1_deg, verbose=True):
    """
    Perform lunar trajectory calculations following the example problem steps
    """
    if verbose:
        print("\n" + "=" * 60)
        print("DETAILED CALCULATION STEPS")
        print("=" * 60)
    
    # Convert angles to radians
    gamma0 = math.radians(gamma0_deg)
    lambda1 = math.radians(lambda1_deg)
    
    # Constants (canonical units based on Earth)
    D = 60.27  # Earth-Moon distance in DU (384,400 km)
    R_s = 10.395  # Moon's sphere of influence radius in DU (66,300 km)
    omega_m = 2.137e-3  # Moon's angular velocity in rad/TU
    miu_earth = 398600.4415  # Earth's gravitational parameter in DU^3/TU^2

    if verbose:
        print(f"Input parameters:")
        print(f"R0 = {R0} DU")
        print(f"V0 = {V0} DU/TU")
        print(f"φ0 (gamma0) = {gamma0_deg}°")
        print(f"λ1 = {lambda1_deg}°")
        print(f"\nConstants:")
        print(f"D = {D} DU")
        print(f"Rs = {R_s} DU")
    
    # Step 1: Calculate specific energy (equation 7.4-2)
    epsilon = (V0**2 / 2) - (1 / R0)
    if verbose:
        print(f"\nStep 1 - Specific Energy (ε):")
        print(f"ε = V0²/2 - 1/R0 = {epsilon:.4f} DU²/TU²")
    
    # Step 2: Calculate specific angular momentum (equation 7.4-3)
    h = R0 * V0 * math.cos(gamma0)
    if verbose:
        print(f"\nStep 2 - Specific Angular Momentum (h):")
        print(f"h = R0 * V0 * cos(γ0) = {h:.4f} DU²/TU")
    
    # Step 3: Calculate r1 at Moon's SOI (equation 7.4-4)
    r1 = math.sqrt(D**2 + R_s**2 - 2*D*R_s*math.cos(lambda1))
    if verbose:
        print(f"\nStep 3 - Distance at Moon's SOI (r1):")
        print(f"r1 = √(D² + Rs² - 2*D*Rs*cos(λ1)) = {r1:.4f} DU")
    
    # Step 4: Calculate v1 at Moon's SOI (equation 7.4-5)
    v1 = math.sqrt(2*(epsilon + h**2/(2*r1**2) + 1/r1))
    if verbose:
        print(f"\nStep 4 - Velocity at Moon's SOI (v1):")
        print(f"v1 = √(2(ε + h²/(2r1²) + 1/r1)) = {v1:.4f} DU/TU")
    
    # Step 5: Calculate φ1 (flight path angle at Moon's SOI) (equation 7.4-6)
    phi1_rad = math.acos(h / (r1 * v1))
    phi1_deg = math.degrees(phi1_rad)
    if verbose:
        print(f"\nStep 5 - Flight Path Angle at Moon's SOI (φ1):")
        print(f"φ1 = arccos(h/(r1*v1)) = {phi1_deg:.2f}°")
    
    # Step 6: Calculate γ1 (phase angle at Moon's SOI) (equation 7.4-7)
    gamma1_rad = math.asin((R_s / r1) * math.sin(lambda1))
    gamma1_deg = math.degrees(gamma1_rad)
    if verbose:
        print(f"\nStep 6 - Phase Angle at Moon's SOI (γ1):")
        print(f"γ1 = arcsin((Rs/r1) * sin(λ1)) = {gamma1_deg:.2f}° = {gamma1_rad:.3f} rad")
    
    # Step 7: Calculate orbital parameters for geocentric trajectory
    if verbose:
        print(f"\nStep 7 - Geocentric Trajectory Parameters:")
    
    # Semi-latus rectum
    p = h**2
    if verbose:
        print(f"p = h² = {p:.3f} DU")
    
    # Semi-major axis
    a = -1 / (2 * epsilon)
    if verbose:
        print(f"a = -1/(2ε) = {a:.2f} DU")
    
    # Eccentricity
    e = math.sqrt(1 + 2*epsilon*h**2)
    if verbose:
        print(f"e = √(1 + 2εh²) = {e:.3f}")
    
    # Step 8: Calculate true anomalies and eccentric anomalies
    if verbose:
        print(f"\nStep 8 - Anomaly Calculations:")
    
    # True anomaly at departure (v0 = 0° since φ0 = 0°)
    v0 = 0
    if verbose:
        print(f"v0 = {v0}° (since φ0 = 0°)")
    
    # True anomaly at arrival
    cos_v1 = (p/r1 - 1) / e
    v1_rad = math.acos(cos_v1)
    v1_deg = math.degrees(v1_rad)
    if verbose:
        print(f"v1 = {v1_deg:.2f}° = {v1_rad:.3f} rad")
    
    # Eccentric anomalies
    E0 = 0  # since v0 = 0°
    E1_rad = 2 * math.atan(math.sqrt((1-e)/(1+e)) * math.tan(v1_rad/2))
    E1_deg = math.degrees(E1_rad)
    if verbose:
        print(f"E0 = {E0}° (since v0 = 0°)")
        print(f"E1 = {E1_deg:.1f}° = {E1_rad:.2f} rad")
    
    # Step 9: Calculate time of flight (TOF)
    if verbose:
        print(f"\nStep 9 - Time of Flight Calculation:")
    sin_E1 = math.sin(E1_rad)
    mean_motion = math.sqrt(1 / a**3)
    tof_tu = (E1_rad - e * sin_E1) / mean_motion
    tof_hours = tof_tu * 0.224114409  # Convert TU to hours (1 TU = 13.44686457 min = 0.224114409 hours)
    
    if verbose:
        print(f"sin E1 = {sin_E1:.3f}")
        print(f"TOF = √(a³) * [(E1 - e*sin E1) - (E0 - e*sin E0)]")
        print(f"TOF = {tof_tu:.2f} TU = {tof_hours:.2f} hours")
    
    # Step 10: Calculate initial phase angle γ0
    if verbose:
        print(f"\nStep 10 - Initial Phase Angle (γ0):")
    gamma0_calc = v1_rad - v0 - gamma1_rad - omega_m * tof_tu
    gamma0_deg_calc = math.degrees(gamma0_calc)
    if verbose:
        print(f"γ0 = v1 - v0 - γ1 - ωm*(t1 - t0)")
        print(f"γ0 = {gamma0_deg_calc:.2f}° = {gamma0_calc:.3f} rad")
    
    return {
        'epsilon': epsilon,
        'h': h,
        'r1': r1,
        'v1': v1,
        'phi1_deg': phi1_deg,
        'gamma1_deg': gamma1_deg,
        'p': p,
        'a': a,
        'e': e,
        'v1_deg': v1_deg,
        'E1_deg': E1_deg,
        'tof_tu': tof_tu,
        'tof_hours': tof_hours,
        'gamma0_calc_deg': gamma0_deg_calc
    }

def lunar_soi_calculations(r1, v1, phi1_deg, lambda1_deg, gamma1_deg, verbose=True):
    """
    Calculate lunar sphere of influence parameters and minimum approach distance
    """
    if verbose:
        print("\n" + "=" * 60)
        print("LUNAR SPHERE OF INFLUENCE CALCULATIONS")
        print("=" * 60)
    
    # Convert to lunar-centered coordinates
    # Constants for Moon
    mu_earth = 398600  # Earth's gravitational parameter (canonical)
    mu_moon = mu_earth / 81.3  # Moon's gravitational parameter
    R_s_km = 66300  # Moon's sphere of influence radius in km
    
    # Convert angles to radians
    phi1_rad = math.radians(phi1_deg)
    lambda1_rad = math.radians(lambda1_deg)
    gamma1_rad = math.radians(gamma1_deg)
    
    # Convert DU/TU to km/s (1 DU = 6378 km, 1 TU = 806.8 min = 13.446 hours)
    v1_kms = v1 * 7.905  # Conversion factor from DU/TU to km/s
    if verbose:
        print(f"v1 = {v1:.4f} DU/TU = {v1_kms:.3f} km/sec")
        print(f"Rs = {10.395} DU = {R_s_km} km")
        print(f"μm = {mu_moon:.6f} (Earth canonical units)")
    
    # Moon's orbital velocity
    vm = 1.018  # km/sec
    if verbose:
        print(f"vm = {vm} km/sec")
    
    # Velocity relative to Moon (equation 7.4-19)
    v2 = math.sqrt(v1_kms**2 + vm**2 - 2*v1_kms*vm*math.cos(phi1_rad))
    if verbose:
        print(f"\nRelative velocity at Moon's SOI:")
        print(f"v2 = √(v1² + vm² - 2*v1*vm*cos(φ1)) = {v2:.3f} km/sec")
    
    # Flight path angle relative to Moon (equation 7.4-20)
    # Calculate the argument for arcsin and check for valid domain
    asin_argument = (vm/v2)*math.cos(lambda1_rad) - (v1_kms/v2)*math.cos(lambda1_rad + gamma1_rad - phi1_rad)
    
    # Check if the argument is within the valid domain [-1, 1] for arcsin
    if not (-1.0 <= asin_argument <= 1.0):
        raise ValueError(f"Invalid flight path geometry for λ₁ = {lambda1_deg}°. "
                        f"arcsin argument = {asin_argument:.6f} is outside [-1, 1] range. "
                        f"This represents a physically impossible trajectory geometry.")
    
    epsilon2_rad = math.asin(asin_argument)
    epsilon2_deg = math.degrees(epsilon2_rad)
    if verbose:
        print(f"ε2 = arcsin[(vm/v2)*cos(λ1) - (v1/v2)*cos(λ1 + γ1 - φ1)] = {epsilon2_deg:.2f}°")
    
    # Lunar-centric trajectory parameters
    if verbose:
        print(f"\nLunar-centric trajectory parameters:")
    
    # Convert μm to km³/sec² for calculations
    mu_m_kms = 4.903e3  # km³/sec²
    
    # Specific angular momentum
    h_lunar = R_s_km * v2 * math.sin(epsilon2_rad)
    if verbose:
        print(f"h = {h_lunar:.0f} km²/sec")
    
    # Specific energy
    epsilon_lunar = (v2**2 / 2) - (mu_m_kms / R_s_km)
    if verbose:
        print(f"ε = v2²/2 - μm/Rs = {epsilon_lunar:.4f}")
    
    # Semi-latus rectum
    p_lunar = (h_lunar**2) / mu_m_kms
    if verbose:
        print(f"p = h²/μm = {p_lunar:.0f} km")
    
    # Eccentricity
    e_lunar = math.sqrt(1 + (2 * epsilon_lunar * h_lunar**2) / mu_m_kms**2)
    if verbose:
        print(f"e = √(1 + 2εh²/μm²) = {e_lunar:.3f}")
    
    # Periapsis distance (minimum approach)
    rp = p_lunar / (1 + e_lunar)
    if verbose:
        print(f"rp = p/(1 + e) = {rp:.0f} km")
    
    # Minimum altitude above lunar surface
    R_moon = 1737  # Moon's radius in km
    hp = rp - R_moon
    if verbose:
        print(f"hp = rp - Rmoon = {hp:.0f} km")
        print(f"\nThis is the minimum distance of approach to the Moon's surface.")
    
    return {
        'v2': v2,
        'epsilon2_deg': epsilon2_deg,
        'h_lunar': h_lunar,
        'epsilon_lunar': epsilon_lunar,
        'p_lunar': p_lunar,
        'e_lunar': e_lunar,
        'rp': rp,
        'hp': hp
    }

def hyperbolic_to_elliptical_conversion(lunar_results, target_perigee_altitude_km=600, verbose=True):
    """
    Convert hyperbolic lunar-centric trajectory to elliptical trajectory
    
    Parameters:
    lunar_results: dictionary with current hyperbolic trajectory parameters
    target_perigee_altitude_km: desired perigee altitude above Moon's surface (default 600 km)
    verbose: whether to print detailed calculations
    
    Returns:
    Dictionary with elliptical trajectory parameters and delta-V requirements
    """
    if verbose:
        print("\n" + "=" * 60)
        print("HYPERBOLIC TO ELLIPTICAL TRAJECTORY CONVERSION")
        print("=" * 60)
    
    # Moon parameters
    R_moon = 1737  # Moon's radius in km
    mu_m_kms = 4.903e3  # Moon's gravitational parameter in km³/sec²
    
    # Extract current hyperbolic trajectory parameters
    rp_hyperbolic = lunar_results['rp']  # perigee radius of hyperbolic trajectory
    e_hyperbolic = lunar_results['e_lunar']  # eccentricity of hyperbolic trajectory
    h_hyperbolic = lunar_results['h_lunar']  # angular momentum of hyperbolic trajectory
    v2_infinity = lunar_results['v2']  # hyperbolic excess velocity
    
    if verbose:
        print(f"Current hyperbolic trajectory:")
        print(f"  Perigee radius = {rp_hyperbolic:.0f} km")
        print(f"  Perigee altitude = {rp_hyperbolic - R_moon:.0f} km above surface")
        print(f"  Eccentricity = {e_hyperbolic:.3f}")
        print(f"  Angular momentum = {h_hyperbolic:.0f} km²/sec")
        print(f"  Hyperbolic excess velocity = {v2_infinity:.3f} km/sec")
    
    # Target elliptical orbit parameters
    rp_elliptical = R_moon + target_perigee_altitude_km  # new perigee radius
    ra_elliptical = rp_hyperbolic  # new apoapsis radius (current perigee becomes apoapsis)
    
    if verbose:
        print(f"\nTarget elliptical trajectory:")
        print(f"  Perigee radius = {rp_elliptical:.0f} km")
        print(f"  Perigee altitude = {target_perigee_altitude_km} km above surface")
        print(f"  Apoapsis radius = {ra_elliptical:.0f} km")
        print(f"  Apoapsis altitude = {ra_elliptical - R_moon:.0f} km above surface")
    
    # Calculate elliptical orbit parameters
    a_elliptical = (rp_elliptical + ra_elliptical) / 2  # semi-major axis
    e_elliptical = (ra_elliptical - rp_elliptical) / (ra_elliptical + rp_elliptical)  # eccentricity
    p_elliptical = a_elliptical * (1 - e_elliptical**2)  # semi-latus rectum
    h_elliptical = math.sqrt(mu_m_kms * p_elliptical)  # new angular momentum
    
    if verbose:
        print(f"\nElliptical orbit parameters:")
        print(f"  Semi-major axis = {a_elliptical:.0f} km")
        print(f"  Eccentricity = {e_elliptical:.3f}")
        print(f"  Semi-latus rectum = {p_elliptical:.0f} km")
        print(f"  Angular momentum = {h_elliptical:.0f} km²/sec")
    
    # Calculate velocities at perigee and apoapsis for elliptical orbit
    v_peri_elliptical = math.sqrt(mu_m_kms * (2/rp_elliptical - 1/a_elliptical))
    v_apo_elliptical = math.sqrt(mu_m_kms * (2/ra_elliptical - 1/a_elliptical))
    
    if verbose:
        print(f"\nElliptical orbit velocities:")
        print(f"  Velocity at perigee = {v_peri_elliptical:.3f} km/sec")
        print(f"  Velocity at apoapsis = {v_apo_elliptical:.3f} km/sec")
    
    # Calculate velocity at perigee of original hyperbolic trajectory using the formula:
    # v_p2 = √(1 + e2) √(μm/r_p2)
    v_peri_hyperbolic = math.sqrt(1 + e_hyperbolic) * math.sqrt(mu_m_kms / rp_hyperbolic)
    
    if verbose:
        print(f"\nOriginal hyperbolic trajectory:")
        print(f"  Velocity at perigee = √(1 + e) √(μm/rp) = {v_peri_hyperbolic:.3f} km/sec")
    
    # Delta-V calculation at the conversion point (current perigee = future apoapsis)
    # The maneuver occurs at the perigee of the hyperbolic trajectory
    delta_v_conversion = v_apo_elliptical - v_peri_hyperbolic
    delta_v_magnitude = abs(delta_v_conversion)
    
    if verbose:
        print(f"\n--- TRAJECTORY CONVERSION MANEUVER ---")
        print(f"Maneuver location: Current perigee (future apoapsis)")
        print(f"Altitude of maneuver: {rp_hyperbolic - R_moon:.0f} km above Moon surface")
        print(f"Delta-V required: {delta_v_magnitude:.3f} km/sec")
        
        if delta_v_conversion < 0:
            print(f"Maneuver type: Retrograde burn (deceleration)")
        else:
            print(f"Maneuver type: Prograde burn (acceleration)")
    
    # Calculate orbital period of elliptical orbit
    T_elliptical = 2 * math.pi * math.sqrt(a_elliptical**3 / mu_m_kms)
    T_hours = T_elliptical / 3600  # convert seconds to hours
    
    if verbose:
        print(f"\nElliptical orbit characteristics:")
        print(f"  Orbital period = {T_elliptical:.0f} seconds = {T_hours:.2f} hours")
    
    # Calculate time from perigee to apoapsis
    t_half_period = T_elliptical / 2
    t_half_hours = t_half_period / 3600
    
    if verbose:
        print(f"  Time from perigee to apoapsis = {t_half_period:.0f} seconds = {t_half_hours:.2f} hours")
    
    # Calculate delta-V for circularization at perigee (600 km altitude)
    if verbose:
        print(f"\n--- CIRCULARIZATION AT PERIGEE (600 km) ---")
    
    # Circular velocity at 600 km altitude
    v_circular_600km = math.sqrt(mu_m_kms / rp_elliptical)
    
    # Delta-V for circularization (difference between circular and elliptical velocities at perigee)
    delta_v_circularization = v_circular_600km - v_peri_elliptical
    delta_v_circ_magnitude = abs(delta_v_circularization)
    
    if verbose:
        print(f"Circular velocity at 600 km altitude = {v_circular_600km:.3f} km/sec")
        print(f"Elliptical velocity at perigee = {v_peri_elliptical:.3f} km/sec")
        print(f"Delta-V for circularization = {delta_v_circ_magnitude:.3f} km/sec")
        
        if delta_v_circularization < 0:
            print(f"Maneuver type: Retrograde burn (deceleration)")
        else:
            print(f"Maneuver type: Prograde burn (acceleration)")
    
    # Total delta-V for complete mission (hyperbolic to circular)
    total_delta_v = delta_v_magnitude + delta_v_circ_magnitude
    if verbose:
        print(f"\n--- TOTAL MISSION DELTA-V ---")
        print(f"Delta-V for hyperbolic to elliptical = {delta_v_magnitude:.3f} km/sec")
        print(f"Delta-V for elliptical to circular = {delta_v_circ_magnitude:.3f} km/sec")
        print(f"TOTAL DELTA-V REQUIRED = {total_delta_v:.3f} km/sec")
    
    return {
        'rp_elliptical': rp_elliptical,
        'ra_elliptical': ra_elliptical,
        'a_elliptical': a_elliptical,
        'e_elliptical': e_elliptical,
        'p_elliptical': p_elliptical,
        'h_elliptical': h_elliptical,
        'v_peri_elliptical': v_peri_elliptical,
        'v_apo_elliptical': v_apo_elliptical,
        'v_peri_hyperbolic': v_peri_hyperbolic,
        'delta_v_conversion': delta_v_conversion,
        'delta_v_magnitude': delta_v_magnitude,
        'orbital_period_hours': T_hours,
        'maneuver_altitude_km': rp_hyperbolic - R_moon,
        'v_circular_600km': v_circular_600km,
        'delta_v_circularization': delta_v_circularization,
        'delta_v_circ_magnitude': delta_v_circ_magnitude,
        'total_delta_v': total_delta_v
    }

def display_organized_mission_summary(inputs, departure_results, geo_results, lunar_results, soi_transit_results, elliptical_results=None):
    """
    Display an organized step-by-step mission summary without redundancies
    """
    print("\n" + "=" * 80)
    print("COMPLETE LUNAR TRANSFER MISSION ANALYSIS")
    print("=" * 80)
    
    # Mission Parameters
    print("\n📋 MISSION PARAMETERS:")
    print(f"   Parking Orbit = Transfer Perigee: {departure_results['altitude_km']:.0f} km altitude")
    print(f"   Transfer Velocity (V0): {inputs['V0']} DU/TU")
    print(f"   Flight Path Angle (γ0): {inputs['gamma0']}°")
    print(f"   Lunar SOI Crossing (λ1): {inputs['lambda1']}°")
    
    # Step 0: Earth Departure
    print("\n🚀 STEP 0: EARTH DEPARTURE")
    print(f"   Parking Orbit Velocity: {departure_results['v_circular_kms']:.3f} km/s")
    print(f"   Transfer Velocity at Departure: {departure_results['v_transfer_kms']:.3f} km/s")
    print(f"   Departure Delta-V: {departure_results['delta_v_departure_kms']:.3f} km/s")
    print(f"   Maneuver Location: {departure_results['altitude_km']:.0f} km altitude")
    
    # Step 1: Earth Departure Trajectory
    print("\n�️ STEP 1: EARTH DEPARTURE TRAJECTORY")
    print(f"   Trajectory Type: {'Elliptical' if geo_results['e'] < 1 else 'Hyperbolic'} (e = {geo_results['e']:.3f})")
    print(f"   Specific Energy: {geo_results['epsilon']:.4f} DU²/TU²")
    print(f"   Angular Momentum: {geo_results['h']:.3f} DU²/TU")
    if geo_results['e'] < 1:
        print(f"   Semi-major Axis: {geo_results['a']:.2f} DU")
    
    # Step 2: Transit to Moon
    print("\n🌍→🌙 STEP 2: EARTH-MOON TRANSIT")
    print(f"   Time of Flight: {geo_results['tof_hours']:.1f} hours ({geo_results['tof_tu']:.2f} TU)")
    print(f"   Distance at Moon's SOI: {geo_results['r1']:.2f} DU")
    print(f"   Velocity at Moon's SOI: {geo_results['v1']:.4f} DU/TU")
    
    # Step 3: Lunar SOI Transit
    print("\n🌙 STEP 3: LUNAR SPHERE OF INFLUENCE APPROACH")
    print(f"   Relative Velocity to Moon: {lunar_results['v2']:.3f} km/s")
    print(f"   Flight Path Angle: {lunar_results['epsilon2_deg']:.1f}°")
    print(f"   Trajectory Type: Hyperbolic (e = {lunar_results['e_lunar']:.3f})")
    print(f"   SOI Entry to Perigee Time: {soi_transit_results['soi_transit_time_hours']:.2f} hours")
    print(f"   Natural Flyby Altitude: {lunar_results['hp']:.0f} km above surface")
    
    if elliptical_results:
        # Step 4: Orbit Insertion
        print("\n🛰️ STEP 4: LUNAR ORBIT INSERTION")
        print(f"   Maneuver Location: {elliptical_results['maneuver_altitude_km']:.0f} km altitude")
        print(f"   Delta-V Required: {elliptical_results['delta_v_magnitude']:.3f} km/s")
        print(f"   Resulting Orbit: {elliptical_results['maneuver_altitude_km']:.0f} × 600 km")
        print(f"   Orbital Period: {elliptical_results['orbital_period_hours']:.2f} hours")
        
        # Step 5: Circularization
        print("\n🔄 STEP 5: ORBIT CIRCULARIZATION (Optional)")
        print(f"   Circularization at: 600 km altitude")
        print(f"   Delta-V Required: {elliptical_results['delta_v_circ_magnitude']:.3f} km/s")
        
        # Mission Summary
        total_mission_time = geo_results['tof_hours'] + soi_transit_results['soi_transit_time_hours']
        total_delta_v_complete = departure_results['delta_v_departure_kms'] + elliptical_results['total_delta_v']
        
        print("\n📊 COMPLETE MISSION SUMMARY:")
        print(f"   Earth Departure Delta-V: {departure_results['delta_v_departure_kms']:.3f} km/s")
        print(f"   Lunar Operations Delta-V: {elliptical_results['total_delta_v']:.3f} km/s")
        print(f"   TOTAL MISSION DELTA-V: {total_delta_v_complete:.3f} km/s")
        print(f"   Earth-Moon Transit Time: {geo_results['tof_hours']:.1f} hours")
        print(f"   SOI Transit Time: {soi_transit_results['soi_transit_time_hours']:.2f} hours")
        print(f"   TOTAL MISSION DURATION: {total_mission_time:.1f} hours")
        print(f"   Final Orbit: 600 km circular")
    
    print("\n" + "=" * 80)

def parametric_study_lambda1():
    """
    Perform a parametric study of lambda1 from 0 to 360 degrees with 5-degree steps
    and plot the relationship between lambda1, total mission delta-V, and time of flight
    """
    print("PARAMETRIC STUDY: LAMBDA1 vs MISSION PARAMETERS")
    print("=" * 60)
    
    # Fixed initial conditions
    R0 = 1.05  # DU - Both parking orbit radius AND transfer trajectory perigee
    V0 = 1.372  # DU/TU
    gamma0 = 0  # degrees
    
    # Lambda1 range: 0 to 360 degrees with 5-degree steps
    lambda1_values = np.arange(0, 361, 5)
    
    # Arrays to store results
    total_delta_v_values = []  # Complete mission delta-V including Earth departure
    time_of_flight_values = []  # Complete mission time including SOI transit
    valid_lambda1_values = []
    
    print(f"Running calculations for {len(lambda1_values)} lambda1 values...")
    print("This may take a moment...")
    
    # Calculate Earth departure delta-V once (it's constant for all lambda1)
    departure_results = calculate_earth_departure_delta_v(R0, V0, verbose=False)
    earth_departure_dv = departure_results['delta_v_departure_kms']
    
    for i, lambda1 in enumerate(lambda1_values):
        try:
            # Perform calculations without verbose output
            geo_results = lunar_trajectory_calculations(R0, V0, gamma0, lambda1, verbose=False)
            lunar_results = lunar_soi_calculations(
                geo_results['r1'], 
                geo_results['v1'], 
                geo_results['phi1_deg'],
                lambda1,
                geo_results['gamma1_deg'],
                verbose=False
            )
            soi_transit_results = calculate_lunar_soi_transit_time(lunar_results, verbose=False)
            elliptical_results = hyperbolic_to_elliptical_conversion(lunar_results, target_perigee_altitude_km=600, verbose=False)
            
            # Store complete mission results
            complete_delta_v = earth_departure_dv + elliptical_results['total_delta_v']
            complete_time = geo_results['tof_hours'] + soi_transit_results['soi_transit_time_hours']
            
            total_delta_v_values.append(complete_delta_v)
            time_of_flight_values.append(complete_time)
            valid_lambda1_values.append(lambda1)
            
            # Progress indicator
            if (i + 1) % 10 == 0:
                print(f"Completed {i + 1}/{len(lambda1_values)} calculations...")
                
        except ValueError as e:
            if "Invalid flight path geometry" in str(e):
                print(f"Skipping λ₁ = {lambda1}° (physically impossible geometry)")
            else:
                print(f"Warning: Calculation failed for lambda1 = {lambda1}°: {str(e)}")
            continue
        except Exception as e:
            print(f"Warning: Unexpected error for lambda1 = {lambda1}°: {str(e)}")
            continue
    
    print(f"Successfully calculated {len(valid_lambda1_values)} data points.")
    skipped_count = len(lambda1_values) - len(valid_lambda1_values)
    if skipped_count > 0:
        print(f"Skipped {skipped_count} lambda1 values due to physically impossible geometries.")
    print(f"Earth departure delta-V (constant): {earth_departure_dv:.3f} km/s")
    
    # Create plots
    create_parametric_plots(valid_lambda1_values, total_delta_v_values, time_of_flight_values)
    
    # Find optimal values
    find_optimal_lambda1(valid_lambda1_values, total_delta_v_values, time_of_flight_values)
    
    return valid_lambda1_values, total_delta_v_values, time_of_flight_values

def create_parametric_plots(lambda1_values, delta_v_values, tof_values):
    """
    Create plots showing the relationship between lambda1 and complete mission parameters
    """
    # Create figure with subplots
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 15))
    
    # Plot 1: Lambda1 vs Complete Mission Delta-V
    ax1.plot(lambda1_values, delta_v_values, 'b-', linewidth=2, marker='o', markersize=3)
    ax1.set_xlabel('Lambda1 (degrees)')
    ax1.set_ylabel('Complete Mission Delta-V (km/s)')
    ax1.set_title('Complete Mission Delta-V vs Lambda1')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 360)
    
    # Add minimum delta-V annotation
    min_delta_v_idx = np.argmin(delta_v_values)
    min_delta_v = delta_v_values[min_delta_v_idx]
    min_lambda1 = lambda1_values[min_delta_v_idx]
    ax1.annotate(f'Min ΔV: {min_delta_v:.3f} km/s\nat λ₁ = {min_lambda1}°', 
                xy=(min_lambda1, min_delta_v), xytext=(min_lambda1 + 50, min_delta_v + 0.5),
                arrowprops=dict(arrowstyle='->', color='red'),
                bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.7))
    
    # Plot 2: Lambda1 vs Complete Mission Time
    ax2.plot(lambda1_values, tof_values, 'r-', linewidth=2, marker='s', markersize=3)
    ax2.set_xlabel('Lambda1 (degrees)')
    ax2.set_ylabel('Complete Mission Time (hours)')
    ax2.set_title('Complete Mission Time vs Lambda1')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, 360)
    
    # Add minimum time annotation
    min_tof_idx = np.argmin(tof_values)
    min_tof = tof_values[min_tof_idx]
    min_tof_lambda1 = lambda1_values[min_tof_idx]
    ax2.annotate(f'Min Time: {min_tof:.1f} hrs\nat λ₁ = {min_tof_lambda1}°', 
                xy=(min_tof_lambda1, min_tof), xytext=(min_tof_lambda1 + 50, min_tof + 20),
                arrowprops=dict(arrowstyle='->', color='red'),
                bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.7))
    
    # Plot 3: Delta-V vs Time of Flight (scatter plot)
    scatter = ax3.scatter(tof_values, delta_v_values, c=lambda1_values, cmap='viridis', 
                         s=50, alpha=0.7, edgecolors='black', linewidth=0.5)
    ax3.set_xlabel('Complete Mission Time (hours)')
    ax3.set_ylabel('Complete Mission Delta-V (km/s)')
    ax3.set_title('Mission Trade-off: Delta-V vs Time of Flight')
    ax3.grid(True, alpha=0.3)
    
    # Add colorbar for lambda1 values
    cbar = plt.colorbar(scatter, ax=ax3)
    cbar.set_label('Lambda1 (degrees)')
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig('lambda1_parametric_study.png', dpi=300, bbox_inches='tight')
    print("Plot saved as 'lambda1_parametric_study.png'")
    plt.show()

def create_complete_parametric_plots(lambda1_values, lunar_delta_v_values, complete_delta_v_values, 
                                   earth_moon_time_values, complete_time_values, earth_departure_dv):
    """
    Create plots showing the complete mission analysis including Earth departure
    """
    # Create figure with subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # Plot 1: Lambda1 vs Complete Mission Delta-V
    ax1.plot(lambda1_values, complete_delta_v_values, 'b-', linewidth=2, marker='o', markersize=3)
    ax1.axhline(y=earth_departure_dv, color='r', linestyle='--', alpha=0.7, 
                label=f'Earth Departure: {earth_departure_dv:.3f} km/s')
    ax1.set_xlabel('Lambda1 (degrees)')
    ax1.set_ylabel('Complete Mission Delta-V (km/s)')
    ax1.set_title('Complete Mission Delta-V vs Lambda1')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 360)
    ax1.legend()
    
    # Add minimum delta-V annotation
    min_delta_v_idx = np.argmin(complete_delta_v_values)
    min_delta_v = complete_delta_v_values[min_delta_v_idx]
    min_lambda1 = lambda1_values[min_delta_v_idx]
    ax1.annotate(f'Min ΔV: {min_delta_v:.3f} km/s\nat λ₁ = {min_lambda1}°', 
                xy=(min_lambda1, min_delta_v), xytext=(min_lambda1 + 50, min_delta_v + 0.5),
                arrowprops=dict(arrowstyle='->', color='red'),
                bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.7))
    
    # Plot 2: Lambda1 vs Complete Mission Time
    ax2.plot(lambda1_values, complete_time_values, 'g-', linewidth=2, marker='s', markersize=3)
    ax2.set_xlabel('Lambda1 (degrees)')
    ax2.set_ylabel('Complete Mission Time (hours)')
    ax2.set_title('Complete Mission Time vs Lambda1')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, 360)
    
    # Add minimum time annotation
    min_time_idx = np.argmin(complete_time_values)
    min_time = complete_time_values[min_time_idx]
    min_time_lambda1 = lambda1_values[min_time_idx]
    ax2.annotate(f'Min Time: {min_time:.1f} hrs\nat λ₁ = {min_time_lambda1}°', 
                xy=(min_time_lambda1, min_time), xytext=(min_time_lambda1 + 50, min_time + 20),
                arrowprops=dict(arrowstyle='->', color='red'),
                bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.7))
    
    # Plot 3: Lunar-only Delta-V vs Lambda1
    ax3.plot(lambda1_values, lunar_delta_v_values, 'm-', linewidth=2, marker='^', markersize=3)
    ax3.set_xlabel('Lambda1 (degrees)')
    ax3.set_ylabel('Lunar Operations Delta-V (km/s)')
    ax3.set_title('Lunar Operations Delta-V vs Lambda1')
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim(0, 360)
    
    # Plot 4: Delta-V vs Time Trade-off (scatter plot)
    scatter = ax4.scatter(complete_time_values, complete_delta_v_values, c=lambda1_values, 
                         cmap='viridis', s=50, alpha=0.7, edgecolors='black', linewidth=0.5)
    ax4.set_xlabel('Complete Mission Time (hours)')
    ax4.set_ylabel('Complete Mission Delta-V (km/s)')
    ax4.set_title('Mission Trade-off: Delta-V vs Time')
    ax4.grid(True, alpha=0.3)
    
    # Add colorbar for lambda1 values
    cbar = plt.colorbar(scatter, ax=ax4)
    cbar.set_label('Lambda1 (degrees)')
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig('complete_mission_parametric_study.png', dpi=300, bbox_inches='tight')
    print("Complete mission plot saved as 'complete_mission_parametric_study.png'")
    plt.show()

def find_optimal_lambda1_complete(lambda1_values, complete_delta_v_values, complete_time_values):
    """
    Find and display optimal lambda1 values for complete mission criteria
    """
    print("\n" + "=" * 60)
    print("OPTIMAL LAMBDA1 ANALYSIS (COMPLETE MISSION)")
    print("=" * 60)
    
    # Convert to numpy arrays for easier manipulation
    lambda1_array = np.array(lambda1_values)
    delta_v_array = np.array(complete_delta_v_values)
    time_array = np.array(complete_time_values)
    
    # Find minimum delta-V
    min_delta_v_idx = np.argmin(delta_v_array)
    min_delta_v = delta_v_array[min_delta_v_idx]
    optimal_lambda1_delta_v = lambda1_array[min_delta_v_idx]
    corresponding_time = time_array[min_delta_v_idx]
    
    print(f"MINIMUM COMPLETE MISSION DELTA-V CRITERION:")
    print(f"  Optimal λ₁ = {optimal_lambda1_delta_v}°")
    print(f"  Minimum ΔV = {min_delta_v:.3f} km/s")
    print(f"  Corresponding Mission Time = {corresponding_time:.1f} hours")
    
    # Find minimum time of flight
    min_time_idx = np.argmin(time_array)
    min_time = time_array[min_time_idx]
    optimal_lambda1_time = lambda1_array[min_time_idx]
    corresponding_delta_v = delta_v_array[min_time_idx]
    
    print(f"\nMINIMUM COMPLETE MISSION TIME CRITERION:")
    print(f"  Optimal λ₁ = {optimal_lambda1_time}°")
    print(f"  Minimum Mission Time = {min_time:.1f} hours")
    print(f"  Corresponding ΔV = {corresponding_delta_v:.3f} km/s")
    
    # Find balanced solution (normalized weighted sum)
    # Normalize both criteria to 0-1 range
    delta_v_normalized = (delta_v_array - np.min(delta_v_array)) / (np.max(delta_v_array) - np.min(delta_v_array))
    time_normalized = (time_array - np.min(time_array)) / (np.max(time_array) - np.min(time_array))
    
    # Equal weighting for balanced solution
    combined_score = 0.5 * delta_v_normalized + 0.5 * time_normalized
    balanced_idx = np.argmin(combined_score)
    balanced_lambda1 = lambda1_array[balanced_idx]
    balanced_delta_v = delta_v_array[balanced_idx]
    balanced_time = time_array[balanced_idx]
    
    print(f"\nBALANCED SOLUTION (Equal weight ΔV and Time):")
    print(f"  Optimal λ₁ = {balanced_lambda1}°")
    print(f"  ΔV = {balanced_delta_v:.3f} km/s")
    print(f"  Mission Time = {balanced_time:.1f} hours")
    
    # Statistics
    print(f"\nSTATISTICS:")
    print(f"  ΔV Range: {np.min(delta_v_array):.3f} - {np.max(delta_v_array):.3f} km/s")
    print(f"  Time Range: {np.min(time_array):.1f} - {np.max(time_array):.1f} hours")
    print(f"  ΔV Standard Deviation: {np.std(delta_v_array):.3f} km/s")
    print(f"  Time Standard Deviation: {np.std(time_array):.1f} hours")
    
    return {
        'min_delta_v': {'lambda1': optimal_lambda1_delta_v, 'delta_v': min_delta_v, 'time': corresponding_time},
        'min_time': {'lambda1': optimal_lambda1_time, 'delta_v': corresponding_delta_v, 'time': min_time},
        'balanced': {'lambda1': balanced_lambda1, 'delta_v': balanced_delta_v, 'time': balanced_time}
    }

def find_optimal_lambda1(lambda1_values, delta_v_values, tof_values):
    """
    Find and display optimal lambda1 values for different criteria
    """
    print("\n" + "=" * 60)
    print("OPTIMAL LAMBDA1 ANALYSIS")
    print("=" * 60)
    
    # Convert to numpy arrays for easier manipulation
    lambda1_array = np.array(lambda1_values)
    delta_v_array = np.array(delta_v_values)
    tof_array = np.array(tof_values)
    
    # Find minimum delta-V
    min_delta_v_idx = np.argmin(delta_v_array)
    min_delta_v = delta_v_array[min_delta_v_idx]
    optimal_lambda1_delta_v = lambda1_array[min_delta_v_idx]
    corresponding_tof = tof_array[min_delta_v_idx]
    
    print(f"MINIMUM DELTA-V CRITERION:")
    print(f"  Optimal λ₁ = {optimal_lambda1_delta_v}°")
    print(f"  Minimum ΔV = {min_delta_v:.3f} km/s")
    print(f"  Corresponding TOF = {corresponding_tof:.1f} hours")
    
    # Find minimum time of flight
    min_tof_idx = np.argmin(tof_array)
    min_tof = tof_array[min_tof_idx]
    optimal_lambda1_tof = lambda1_array[min_tof_idx]
    corresponding_delta_v = delta_v_array[min_tof_idx]
    
    print(f"\nMINIMUM TIME OF FLIGHT CRITERION:")
    print(f"  Optimal λ₁ = {optimal_lambda1_tof}°")
    print(f"  Minimum TOF = {min_tof:.1f} hours")
    print(f"  Corresponding ΔV = {corresponding_delta_v:.3f} km/s")
    
    # Find balanced solution (normalized weighted sum)
    # Normalize both criteria to 0-1 range
    delta_v_normalized = (delta_v_array - np.min(delta_v_array)) / (np.max(delta_v_array) - np.min(delta_v_array))
    tof_normalized = (tof_array - np.min(tof_array)) / (np.max(tof_array) - np.min(tof_array))
    
    # Equal weighting for balanced solution
    combined_score = 0.5 * delta_v_normalized + 0.5 * tof_normalized
    balanced_idx = np.argmin(combined_score)
    balanced_lambda1 = lambda1_array[balanced_idx]
    balanced_delta_v = delta_v_array[balanced_idx]
    balanced_tof = tof_array[balanced_idx]
    
    print(f"\nBALANCED SOLUTION (Equal weight ΔV and TOF):")
    print(f"  Optimal λ₁ = {balanced_lambda1}°")
    print(f"  ΔV = {balanced_delta_v:.3f} km/s")
    print(f"  TOF = {balanced_tof:.1f} hours")
    
    # Statistics
    print(f"\nSTATISTICS:")
    print(f"  ΔV Range: {np.min(delta_v_array):.3f} - {np.max(delta_v_array):.3f} km/s")
    print(f"  TOF Range: {np.min(tof_array):.1f} - {np.max(tof_array):.1f} hours")
    print(f"  ΔV Standard Deviation: {np.std(delta_v_array):.3f} km/s")
    print(f"  TOF Standard Deviation: {np.std(tof_array):.1f} hours")
    
    return {
        'min_delta_v': {'lambda1': optimal_lambda1_delta_v, 'delta_v': min_delta_v, 'tof': corresponding_tof},
        'min_tof': {'lambda1': optimal_lambda1_tof, 'delta_v': corresponding_delta_v, 'tof': min_tof},
        'balanced': {'lambda1': balanced_lambda1, 'delta_v': balanced_delta_v, 'tof': balanced_tof}
    }

def main():
    """
    Main function to run the complete lunar transfer trajectory analysis
    """
    print("LUNAR TRANSFER TRAJECTORY ANALYSIS PROGRAM")
    print("=" * 60)
    
    # Ask user what analysis to perform
    print("Choose analysis type:")
    print("1. Single calculation with user-defined lambda1")
    print("2. Parametric study: lambda1 from 0-360° (5° steps)")
    
    while True:
        try:
            choice = input("Enter choice (1 or 2): ").strip()
            if choice in ['1', '2']:
                break
            else:
                print("Please enter 1 or 2")
        except:
            print("Please enter 1 or 2")
    
    if choice == '1':
        # Original single calculation
        user_inputs = get_user_input()
        
        # Extract variables
        R0 = user_inputs['R0']
        V0 = user_inputs['V0']
        gamma0 = user_inputs['gamma0']
        lambda1 = user_inputs['lambda1']
        
        # Perform complete mission analysis with detailed output
        departure_results = calculate_earth_departure_delta_v(R0, V0, verbose=True)
        geo_results = lunar_trajectory_calculations(R0, V0, gamma0, lambda1, verbose=True)
        lunar_results = lunar_soi_calculations(
            geo_results['r1'], 
            geo_results['v1'], 
            geo_results['phi1_deg'],
            lambda1,
            geo_results['gamma1_deg'],
            verbose=True
        )
        soi_transit_results = calculate_lunar_soi_transit_time(lunar_results, verbose=True)
        elliptical_results = hyperbolic_to_elliptical_conversion(lunar_results, target_perigee_altitude_km=600, verbose=True)
        
        # Display complete organized mission summary
        display_organized_mission_summary(user_inputs, departure_results, geo_results, lunar_results, 
                                        soi_transit_results, elliptical_results)
        
        return user_inputs, departure_results, geo_results, lunar_results, soi_transit_results, elliptical_results
    
    else:
        # Parametric study
        lambda1_values, delta_v_values, tof_values = parametric_study_lambda1()
        optimal_results = find_optimal_lambda1(lambda1_values, delta_v_values, tof_values)
        
        return lambda1_values, delta_v_values, tof_values, optimal_results

if __name__ == "__main__":
    results = main()
