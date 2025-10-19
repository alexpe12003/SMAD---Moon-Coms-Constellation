
"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This function computes the classical orbital elements (coe)
from the state vector (R,V) using Algorithm 4.1.

  mu   - gravitational parameter (km^3/s^2)
  R    - position vector in the geocentric equatorial frame (km)
  V    - velocity vector in the geocentric equatorial frame (km)
  r, v - the magnitudes of R and V
  vr   - radial velocity component (km/s)
  H    - the angular momentum vector (km^2/s)
  h    - the magnitude of H (km^2/s)
  incl - inclination of the orbit (rad)
  N    - the node line vector (km^2/s)
  n    - the magnitude of N
  cp   - cross product of N and R
  RA   - right ascension of the ascending node (rad)
  E    - eccentricity vector
  e    - eccentricity (magnitude of E)
  eps  - a small number below which the eccentricity is considered
         to be zero
  w    - argument of perigee (rad)
  TA   - true anomaly (rad)
  a    - semimajor axis (km)
  pi   - 3.1415926...
  coe  - vector of orbital elements [h e RA incl w TA a]
    Note: All angles are in radians
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

import numpy as np

def coe_from_sv(R, V, mu):
    """
    Compute classical orbital elements from state vector (R, V).
    Returns: [h, e, RA, incl, w, TA, a]
    """
    eps = 1e-6  # Small threshold for zero eccentricity

    r = np.linalg.norm(R)
    v = np.linalg.norm(V)

    vr = np.dot(R, V) / r  # Radial velocity component

    H = np.cross(R, V)     # Angular momentum vector
    h = np.linalg.norm(H)

    # Inclination (rad)
    incl = np.arccos(H[2] / h)

    # Node line vector
    K = np.array([0, 0, 1])
    N = np.cross(K, H)
    n = np.linalg.norm(N)

    # Right ascension of ascending node (RA)
    if not np.isclose(incl, 0):  # Inclined orbit
        RA = np.arccos(N[0] / n)
        if N[1] < 0:
            RA = 2 * np.pi - RA
    else:  # Equatorial orbit
        RA = 0.0

    # Eccentricity vector and magnitude
    E = (1 / mu) * ((v ** 2 - mu / r) * R - r * vr * V)
    e = np.linalg.norm(E)

    # Argument of perigee (w)
    if not np.isclose(incl, 0):  # Inclined orbit
        if e > eps:  # Non-circular
            w = np.arccos(np.dot(N, E) / (n * e))
            if E[2] < 0:
                w = 2 * np.pi - w
        else:  # Circular
            w = 0.0
    else:  # Equatorial orbit
        if e > eps:  # Non-circular
            w = np.arccos(E[0] / e)
            if E[1] < 0:
                w = 2 * np.pi - w
        else:  # Circular
            w = 0.0

    # True anomaly (TA)
    if not np.isclose(incl, 0):  # Inclined orbit
        if e > eps:  # Non-circular
            TA = np.arccos(np.dot(E, R) / (e * r))
            if vr < 0:
                TA = 2 * np.pi - TA
        else:  # Circular
            TA = np.arccos(np.dot(N, R) / (n * r))
            if R[2] < 0:
                TA = 2 * np.pi - TA
    else:  # Equatorial orbit
        if e > eps:  # Non-circular
            TA = np.arccos(np.dot(E, R) / (e * r))
            if vr < 0:
                TA = 2 * np.pi - TA
        else:  # Circular
            TA = np.arccos(R[0] / r)
            if R[1] < 0:
                TA = 2 * np.pi - TA

    # Semimajor axis (a)
    # For elliptical orbits, a > 0; for hyperbolic, a < 0
    a = h ** 2 / mu / (1 - e ** 2)

    # Return orbital elements as a list
    # h: specific angular momentum
    # e: eccentricity
    # RA: right ascension of ascending node
    # incl: inclination
    # w: argument of perigee
    # TA: true anomaly
    # a: semimajor axis
    coe = [h, e, RA, incl, w, TA, a]
    return coe
    # Step-by-step explanation:
    # 1. Compute magnitudes of position and velocity vectors
    # 2. Compute radial velocity (component of velocity along position vector)
    # 3. Compute angular momentum vector and its magnitude
    # 4. Compute inclination (angle between angular momentum and z-axis)
    # 5. Compute node line (intersection of orbital plane and equatorial plane)
    # 6. Compute right ascension of ascending node (angle from x-axis to node line)
    # 7. Compute eccentricity vector and its magnitude
    # 8. Compute argument of perigee (angle from node line to eccentricity vector)
    # 9. Compute true anomaly (angle from eccentricity vector to position vector)
    # 10. Compute semimajor axis (size of the orbit)

if __name__ == "__main__":
    # Example test case: Earth orbit
    # Position and velocity vectors in km and km/s
    R = np.array([7000.0, -12124.0, 0.0])
    V = np.array([2.6679, 4.6210, 0.0])
    mu = 398600.4418  # Earth's gravitational parameter (km^3/s^2)

    coe = coe_from_sv(R, V, mu)
    print("State Variables:")
    print(f"R (km): {R}")
    print(f"V (km/s): {V}")
    print("Classical Orbital Elements:")
    print(f"h (km^2/s): {coe[0]}")
    print(f"e (eccentricity): {coe[1]}")
    print(f"RA (rad): {coe[2]}")
    print(f"incl (rad): {coe[3]}")
    print(f"w (rad): {coe[4]}")
    print(f"TA (rad): {coe[5]}")
    print(f"a (km): {coe[6]}")

