
"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This function computes the state vector (r,v) from the
classical orbital elements (coe).

  mu   - gravitational parameter (km^3;s^2)
  coe  - orbital elements [h e RA incl w TA]
         where
             h    = angular momentum (km^2/s)
             e    = eccentricity
             RA   = right ascension of the ascending node (rad)
             incl = inclination of the orbit (rad)
             w    = argument of perigee (rad)
             TA   = true anomaly (rad)
  R3_w - Rotation matrix about the z-axis through the angle w
  R1_i - Rotation matrix about the x-axis through the angle i
  R3_W - Rotation matrix about the z-axis through the angle RA
  Q_pX - Matrix of the transformation from perifocal to geocentric 
         equatorial frame
  rp   - position vector in the perifocal frame (km)
  vp   - velocity vector in the perifocal frame (km/s)
  r    - position vector in the geocentric equatorial frame (km)
  v    - velocity vector in the geocentric equatorial frame (km/s)

  User M-functions required: none
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

import numpy as np

def sv_from_coe(coe, mu):
    """
    Compute state vector (r, v) from classical orbital elements.
    coe: [h, e, RA, incl, w, TA]
    mu: gravitational parameter (km^3/s^2)
    Returns: r, v (row vectors in geocentric equatorial frame)
    """
    # Unpack classical orbital elements
    h, e, RA, incl, w, TA = coe

    # Compute position vector in perifocal frame (rp)
    # rp = (h^2/mu) * (1/(1 + e*cos(TA))) * [cos(TA), sin(TA), 0]
    rp = (h**2 / mu) * (1 / (1 + e * np.cos(TA))) * (
        np.cos(TA) * np.array([1, 0, 0]) + np.sin(TA) * np.array([0, 1, 0])
    )

    # Compute velocity vector in perifocal frame (vp)
    # vp = (mu/h) * [-sin(TA), e + cos(TA), 0]
    vp = (mu / h) * (
        -np.sin(TA) * np.array([1, 0, 0]) + (e + np.cos(TA)) * np.array([0, 1, 0])
    )

    # Rotation matrix about z-axis by RA (right ascension of ascending node)
    R3_W = np.array([
        [np.cos(RA),  np.sin(RA), 0],
        [-np.sin(RA), np.cos(RA), 0],
        [0, 0, 1]
    ])

    # Rotation matrix about x-axis by inclination
    R1_i = np.array([
        [1, 0, 0],
        [0, np.cos(incl), np.sin(incl)],
        [0, -np.sin(incl), np.cos(incl)]
    ])

    # Rotation matrix about z-axis by argument of perigee (w)
    R3_w = np.array([
        [np.cos(w),  np.sin(w), 0],
        [-np.sin(w), np.cos(w), 0],
        [0, 0, 1]
    ])

    # Transformation matrix from perifocal to geocentric equatorial frame
    # Q_pX = (R3_w * R1_i * R3_W)'
    # The .T is the transpose, matching MATLAB's '
    Q_pX = (R3_w @ R1_i @ R3_W).T

    # Transform position and velocity to geocentric equatorial frame
    r = Q_pX @ rp  # Position vector in geocentric equatorial frame
    v = Q_pX @ vp  # Velocity vector in geocentric equatorial frame

    # Return as row vectors (1D numpy arrays)
    return r, v

if __name__ == "__main__":
    # Example test case: Earth orbit
    # Orbital elements: [h, e, RA, incl, w, TA]
    # All angles in radians
    h = 64692.619         # km^2/s, specific angular momentum
    e = 0.49999           # eccentricity
    RA = np.deg2rad(0)    # right ascension of ascending node (radians)
    incl = np.deg2rad(0)  # inclination (radians)
    w = 1.0472492648      # argument of perigee (radians)
    TA = 4.1887411933     # true anomaly (radians)
    mu = 398600.4418      # Earth's gravitational parameter (km^3/s^2)

    coe = [h, e, RA, incl, w, TA]
    r, v = sv_from_coe(coe, mu)
    print("Position vector r (km):", r)
    print("Velocity vector v (km/s):", v)