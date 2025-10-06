import math

def get_user_input():
    """
    Get 4 orbital mechanics variables from user input: 
    R0 (initial position in DU), V0 (initial velocity in DU/TU), 
    gamma0 (initial flight path angle in degrees), and lambda1 (lunar sphere crossing point in degrees)
    Returns a dictionary with the input values
    """
    print("Please enter the following parameters:")
    print("-" * 50)
    
    variables = {}
    
    # Define the variables to collect
    var_names = ['R0', 'V0', 'gamma0', 'lambda1']
    var_descriptions = [
        'R0 (Initial position in DU)',
        'V0 (Initial velocity in DU/TU)', 
        'gamma0 (Initial flight path angle in degrees)',
        'lambda1 (Point where geocentric trajectory crosses lunar sphere of influence in degrees)'
    ]
    
    # Collect each variable with input validation
    for var_name, description in zip(var_names, var_descriptions):
        while True:
            try:
                value = float(input(f"Enter {description}: "))
                variables[var_name] = value
                break
            except ValueError:
                print("Error: Please enter a valid number.")
                continue
    
    return variables

def lunar_trajectory_calculations(R0, V0, gamma0_deg, lambda1_deg):
    """
    Perform lunar trajectory calculations following the example problem steps
    """
    print("\n" + "=" * 60)
    print("LUNAR TRAJECTORY CALCULATIONS")
    print("=" * 60)
    
    # Convert angles to radians
    gamma0 = math.radians(gamma0_deg)
    lambda1 = math.radians(lambda1_deg)
    
    # Constants (canonical units based on Earth)
    D = 60.27  # Earth-Moon distance in DU (384,400 km)
    R_s = 10.395  # Moon's sphere of influence radius in DU (66,300 km)
    omega_m = 2.137e-3  # Moon's angular velocity in rad/TU
    miu_earth = 398600.4415  # Earth's gravitational parameter in DU^3/TU^2

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
    print(f"\nStep 1 - Specific Energy (ε):")
    print(f"ε = V0²/2 - 1/R0 = {epsilon:.4f} DU²/TU²")
    
    # Step 2: Calculate specific angular momentum (equation 7.4-3)
    h = R0 * V0 * math.cos(gamma0)
    print(f"\nStep 2 - Specific Angular Momentum (h):")
    print(f"h = R0 * V0 * cos(γ0) = {h:.4f} DU²/TU")
    
    # Step 3: Calculate r1 at Moon's SOI (equation 7.4-4)
    r1 = math.sqrt(D**2 + R_s**2 - 2*D*R_s*math.cos(lambda1))
    print(f"\nStep 3 - Distance at Moon's SOI (r1):")
    print(f"r1 = √(D² + Rs² - 2*D*Rs*cos(λ1)) = {r1:.4f} DU")
    
    # Step 4: Calculate v1 at Moon's SOI (equation 7.4-5)
    v1 = math.sqrt(2*(epsilon + h**2/(2*r1**2) + 1/r1))
    print(f"\nStep 4 - Velocity at Moon's SOI (v1):")
    print(f"v1 = √(2(ε + h²/(2r1²) + 1/r1)) = {v1:.4f} DU/TU")
    
    # Step 5: Calculate φ1 (flight path angle at Moon's SOI) (equation 7.4-6)
    phi1_rad = math.acos(h / (r1 * v1))
    phi1_deg = math.degrees(phi1_rad)
    print(f"\nStep 5 - Flight Path Angle at Moon's SOI (φ1):")
    print(f"φ1 = arccos(h/(r1*v1)) = {phi1_deg:.2f}°")
    
    # Step 6: Calculate γ1 (phase angle at Moon's SOI) (equation 7.4-7)
    gamma1_rad = math.asin((R_s / r1) * math.sin(lambda1))
    gamma1_deg = math.degrees(gamma1_rad)
    print(f"\nStep 6 - Phase Angle at Moon's SOI (γ1):")
    print(f"γ1 = arcsin((Rs/r1) * sin(λ1)) = {gamma1_deg:.2f}° = {gamma1_rad:.3f} rad")
    
    # Step 7: Calculate orbital parameters for geocentric trajectory
    print(f"\nStep 7 - Geocentric Trajectory Parameters:")
    
    # Semi-latus rectum
    p = h**2
    print(f"p = h² = {p:.3f} DU")
    
    # Semi-major axis
    a = -1 / (2 * epsilon)
    print(f"a = -1/(2ε) = {a:.2f} DU")
    
    # Eccentricity
    e = math.sqrt(1 + 2*epsilon*h**2)
    print(f"e = √(1 + 2εh²) = {e:.3f}")
    
    # Step 8: Calculate true anomalies and eccentric anomalies
    print(f"\nStep 8 - Anomaly Calculations:")
    
    # True anomaly at departure (v0 = 0° since φ0 = 0°)
    v0 = 0
    print(f"v0 = {v0}° (since φ0 = 0°)")
    
    # True anomaly at arrival
    cos_v1 = (p/r1 - 1) / e
    v1_rad = math.acos(cos_v1)
    v1_deg = math.degrees(v1_rad)
    print(f"v1 = {v1_deg:.2f}° = {v1_rad:.3f} rad")
    
    # Eccentric anomalies
    E0 = 0  # since v0 = 0°
    E1_rad = 2 * math.atan(math.sqrt((1-e)/(1+e)) * math.tan(v1_rad/2))
    E1_deg = math.degrees(E1_rad)
    print(f"E0 = {E0}° (since v0 = 0°)")
    print(f"E1 = {E1_deg:.1f}° = {E1_rad:.2f} rad")
    
    # Step 9: Calculate time of flight (TOF)
    print(f"\nStep 9 - Time of Flight Calculation:")
    sin_E1 = math.sin(E1_rad)
    mean_motion = math.sqrt(1 / a**3)
    tof_tu = (E1_rad - e * sin_E1) / mean_motion
    tof_hours = tof_tu * 0.224114409  # Convert TU to hours (1 TU = 13.44686457 min = 0.224114409 hours)
    
    print(f"sin E1 = {sin_E1:.3f}")
    print(f"TOF = √(a³) * [(E1 - e*sin E1) - (E0 - e*sin E0)]")
    print(f"TOF = {tof_tu:.2f} TU = {tof_hours:.2f} hours")
    
    # Step 10: Calculate initial phase angle γ0
    print(f"\nStep 10 - Initial Phase Angle (γ0):")
    gamma0_calc = v1_rad - v0 - gamma1_rad - omega_m * tof_tu
    gamma0_deg_calc = math.degrees(gamma0_calc)
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

def lunar_soi_calculations(r1, v1, phi1_deg, lambda1_deg, gamma1_deg):
    """
    Calculate lunar sphere of influence parameters and minimum approach distance
    """
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
    print(f"v1 = {v1:.4f} DU/TU = {v1_kms:.3f} km/sec")
    print(f"Rs = {10.395} DU = {R_s_km} km")
    print(f"μm = {mu_moon:.6f} (Earth canonical units)")
    
    # Moon's orbital velocity
    vm = 1.018  # km/sec
    print(f"vm = {vm} km/sec")
    
    # Velocity relative to Moon (equation 7.4-19)
    v2 = math.sqrt(v1_kms**2 + vm**2 - 2*v1_kms*vm*math.cos(phi1_rad))
    print(f"\nRelative velocity at Moon's SOI:")
    print(f"v2 = √(v1² + vm² - 2*v1*vm*cos(φ1)) = {v2:.3f} km/sec")
    
    # Flight path angle relative to Moon (equation 7.4-20)
    epsilon2_rad = math.asin((vm/v2)*math.cos(lambda1_rad) - (v1_kms/v2)*math.cos(lambda1_rad + gamma1_rad - phi1_rad))
    epsilon2_deg = math.degrees(epsilon2_rad)
    print(f"ε2 = arcsin[(vm/v2)*cos(λ1) - (v1/v2)*cos(λ1 + γ1 - φ1)] = {epsilon2_deg:.2f}°")
    
    # Lunar-centric trajectory parameters
    print(f"\nLunar-centric trajectory parameters:")
    
    # Convert μm to km³/sec² for calculations
    mu_m_kms = 4.903e3  # km³/sec²
    
    # Specific angular momentum
    h_lunar = R_s_km * v2 * math.sin(epsilon2_rad)
    print(f"h = {h_lunar:.0f} km²/sec")
    
    # Specific energy
    epsilon_lunar = (v2**2 / 2) - (mu_m_kms / R_s_km)
    print(f"ε = v2²/2 - μm/Rs = {epsilon_lunar:.4f}")
    
    # Semi-latus rectum
    p_lunar = (h_lunar**2) / mu_m_kms
    print(f"p = h²/μm = {p_lunar:.0f} km")
    
    # Eccentricity
    e_lunar = math.sqrt(1 + (2 * epsilon_lunar * h_lunar**2) / mu_m_kms**2)
    print(f"e = √(1 + 2εh²/μm²) = {e_lunar:.3f}")
    
    # Periapsis distance (minimum approach)
    rp = p_lunar / (1 + e_lunar)
    print(f"rp = p/(1 + e) = {rp:.0f} km")
    
    # Minimum altitude above lunar surface
    R_moon = 1737  # Moon's radius in km
    hp = rp - R_moon
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

def hyperbolic_to_elliptical_conversion(lunar_results, target_perigee_altitude_km=600):
    """
    Convert hyperbolic lunar-centric trajectory to elliptical trajectory
    
    Parameters:
    lunar_results: dictionary with current hyperbolic trajectory parameters
    target_perigee_altitude_km: desired perigee altitude above Moon's surface (default 600 km)
    
    Returns:
    Dictionary with elliptical trajectory parameters and delta-V requirements
    """
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
    
    print(f"Current hyperbolic trajectory:")
    print(f"  Perigee radius = {rp_hyperbolic:.0f} km")
    print(f"  Perigee altitude = {rp_hyperbolic - R_moon:.0f} km above surface")
    print(f"  Eccentricity = {e_hyperbolic:.3f}")
    print(f"  Angular momentum = {h_hyperbolic:.0f} km²/sec")
    print(f"  Hyperbolic excess velocity = {v2_infinity:.3f} km/sec")
    
    # Target elliptical orbit parameters
    rp_elliptical = R_moon + target_perigee_altitude_km  # new perigee radius
    ra_elliptical = rp_hyperbolic  # new apoapsis radius (current perigee becomes apoapsis)
    
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
    
    print(f"\nElliptical orbit parameters:")
    print(f"  Semi-major axis = {a_elliptical:.0f} km")
    print(f"  Eccentricity = {e_elliptical:.3f}")
    print(f"  Semi-latus rectum = {p_elliptical:.0f} km")
    print(f"  Angular momentum = {h_elliptical:.0f} km²/sec")
    
    # Calculate velocities at perigee and apoapsis for elliptical orbit
    v_peri_elliptical = math.sqrt(mu_m_kms * (2/rp_elliptical - 1/a_elliptical))
    v_apo_elliptical = math.sqrt(mu_m_kms * (2/ra_elliptical - 1/a_elliptical))
    
    print(f"\nElliptical orbit velocities:")
    print(f"  Velocity at perigee = {v_peri_elliptical:.3f} km/sec")
    print(f"  Velocity at apoapsis = {v_apo_elliptical:.3f} km/sec")
    
    # Calculate velocity at perigee of original hyperbolic trajectory using the formula:
    # v_p2 = √(1 + e2) √(μm/r_p2)
    v_peri_hyperbolic = math.sqrt(1 + e_hyperbolic) * math.sqrt(mu_m_kms / rp_hyperbolic)
    
    print(f"\nOriginal hyperbolic trajectory:")
    print(f"  Velocity at perigee = √(1 + e) √(μm/rp) = {v_peri_hyperbolic:.3f} km/sec")
    
    # Delta-V calculation at the conversion point (current perigee = future apoapsis)
    # The maneuver occurs at the perigee of the hyperbolic trajectory
    delta_v_conversion = v_apo_elliptical - v_peri_hyperbolic
    delta_v_magnitude = abs(delta_v_conversion)
    
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
    
    print(f"\nElliptical orbit characteristics:")
    print(f"  Orbital period = {T_elliptical:.0f} seconds = {T_hours:.2f} hours")
    
    # Calculate time from perigee to apoapsis
    t_half_period = T_elliptical / 2
    t_half_hours = t_half_period / 3600
    
    print(f"  Time from perigee to apoapsis = {t_half_period:.0f} seconds = {t_half_hours:.2f} hours")
    
    # Calculate delta-V for circularization at perigee (600 km altitude)
    print(f"\n--- CIRCULARIZATION AT PERIGEE (600 km) ---")
    
    # Circular velocity at 600 km altitude
    v_circular_600km = math.sqrt(mu_m_kms / rp_elliptical)
    
    # Delta-V for circularization (difference between circular and elliptical velocities at perigee)
    delta_v_circularization = v_circular_600km - v_peri_elliptical
    delta_v_circ_magnitude = abs(delta_v_circularization)
    
    print(f"Circular velocity at 600 km altitude = {v_circular_600km:.3f} km/sec")
    print(f"Elliptical velocity at perigee = {v_peri_elliptical:.3f} km/sec")
    print(f"Delta-V for circularization = {delta_v_circ_magnitude:.3f} km/sec")
    
    if delta_v_circularization < 0:
        print(f"Maneuver type: Retrograde burn (deceleration)")
    else:
        print(f"Maneuver type: Prograde burn (acceleration)")
    
    # Total delta-V for complete mission (hyperbolic to circular)
    total_delta_v = delta_v_magnitude + delta_v_circ_magnitude
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

def display_summary(inputs, geo_results, lunar_results, elliptical_results=None):
    """
    Display a comprehensive summary of all calculations
    """
    print("\n" + "=" * 60)
    print("CALCULATION SUMMARY")
    print("=" * 60)
    
    print("INPUT PARAMETERS:")
    print(f"  R0 = {inputs['R0']} DU")
    print(f"  V0 = {inputs['V0']} DU/TU") 
    print(f"  γ0 = {inputs['gamma0']}°")
    print(f"  λ1 = {inputs['lambda1']}°")
    
    print("\nGEOCENTRIC TRAJECTORY RESULTS:")
    print(f"  Specific Energy (ε) = {geo_results['epsilon']:.4f} DU²/TU²")
    print(f"  Angular Momentum (h) = {geo_results['h']:.3f} DU²/TU")
    print(f"  Semi-major axis (a) = {geo_results['a']:.2f} DU")
    print(f"  Eccentricity (e) = {geo_results['e']:.3f}")
    print(f"  Time of Flight = {geo_results['tof_hours']:.1f} hours")
    print(f"  Calculated γ0 = {geo_results['gamma0_calc_deg']:.1f}°")
    
    print("\nLUNAR APPROACH RESULTS:")
    print(f"  Velocity at Moon's SOI = {geo_results['v1']:.4f} DU/TU")
    print(f"  Flight path angle at SOI = {geo_results['phi1_deg']:.1f}°")
    print(f"  Relative velocity to Moon = {lunar_results['v2']:.3f} km/sec")
    print(f"  Minimum approach distance (Perigee Radius) = {lunar_results['hp']:.0f} km above surface")
    
    if elliptical_results:
        print("\nELLIPTICAL ORBIT CONVERSION:")
        print(f"  Target perigee altitude = 600 km")
        print(f"  Target apoapsis altitude = {elliptical_results['maneuver_altitude_km']:.0f} km")
        print(f"  Elliptical orbit eccentricity = {elliptical_results['e_elliptical']:.3f}")
        print(f"  Orbital period = {elliptical_results['orbital_period_hours']:.2f} hours")
        print(f"  Delta-V for conversion = {elliptical_results['delta_v_magnitude']:.3f} km/sec")
        print(f"  Delta-V for circularization = {elliptical_results['delta_v_circ_magnitude']:.3f} km/sec")
        print(f"  TOTAL MISSION DELTA-V = {elliptical_results['total_delta_v']:.3f} km/sec")
    
    print("=" * 60)

def main():
    """
    Main function to run the complete lunar transfer trajectory analysis
    """
    print("LUNAR TRANSFER TRAJECTORY ANALYSIS PROGRAM")
    print("=" * 60)
    
    # Get user inputs
    user_inputs = get_user_input()
    
    # Extract variables
    R0 = user_inputs['R0']
    V0 = user_inputs['V0']
    gamma0 = user_inputs['gamma0']
    lambda1 = user_inputs['lambda1']
    
    # Perform geocentric trajectory calculations
    geo_results = lunar_trajectory_calculations(R0, V0, gamma0, lambda1)
    
    # Perform lunar SOI calculations
    lunar_results = lunar_soi_calculations(
        geo_results['r1'], 
        geo_results['v1'], 
        geo_results['phi1_deg'],
        lambda1,
        geo_results['gamma1_deg']
    )
    
    # Perform hyperbolic to elliptical trajectory conversion
    elliptical_results = hyperbolic_to_elliptical_conversion(lunar_results, target_perigee_altitude_km=600)
    
    # Display comprehensive summary
    display_summary(user_inputs, geo_results, lunar_results, elliptical_results)
    
    return user_inputs, geo_results, lunar_results, elliptical_results

if __name__ == "__main__":
    inputs, geo_results, lunar_results, elliptical_results = main()
