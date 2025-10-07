import numpy as np  # Import NumPy for vector and math operations
# Lambert's problem solver using universal variable formulation and knowing r1 and r2 position Vectors module, time of flight t, and delta theta
def lambert(r1, r2, t, delta_theta, string, mu):
    """
    Solves Lambert's problem for two position vectors and time of flight.

    Parameters:
        R1 (np.ndarray): Initial position vector (km)
        R2 (np.ndarray): Final position vector (km)
        t (float): Time of flight from R1 to R2 (s)
        string (str): 'pro' for prograde, 'retro' for retrograde (default 'pro')
        mu (float): Gravitational parameter (km^3/s^2, default is Earth's)

    Returns:
        V1 (np.ndarray): Initial velocity vector (km/s)
        V2 (np.ndarray): Final velocity vector (km/s)
    """

    # Construct R1 and R2 from r1, r2, theta1, theta2
    theta1 = 0.0  # starting true anomaly (can be parameterized)
    theta2 = delta_theta  # final true anomaly (can be parameterized)
    R1 = np.array([r1 * np.cos(theta1), r1 * np.sin(theta1), 0])
    R2 = np.array([r2 * np.cos(theta2), r2 * np.sin(theta2), 0])
    # Angle between R1 and R2 (radians)
    theta = theta2 - theta1

    # Calculate the cross product and angle between R1 and R2
    c12 = np.cross(R1, R2)  # Used to determine orbit direction
  

    # Determine if the transfer is prograde or retrograde
    # Prograde: motion in the direction of planetary rotation
    # Retrograde: opposite direction
    if string not in ['pro', 'retro']:
        string = 'pro'  # Default to prograde if not specified
        print('\n ** Prograde trajectory assumed.\n')

    # Adjust theta based on orbit direction and cross product sign
    if string == 'pro':
        if c12[2] <= 0:
            theta = 2 * np.pi - theta
    elif string == 'retro':
        if c12[2] >= 0:
            theta = 2 * np.pi - theta

    # Compute the constant A (geometry of the transfer)
    # Equation 5.35 from Bate, Mueller, White
    # A relates the geometry of the transfer orbit to the time of flight
    A = np.sin(theta) * np.sqrt(r1 * r2 / (1 - np.cos(theta)))

    # Universal variable formulation: solve for z using Newton's method
    z = 0.0  # Initial guess for universal anomaly
    tol = 1e-8  # Tolerance for Newton's method convergence
    nmax = 5000  # Maximum number of iterations
    ratio = 1  # Newton's method update ratio
    n = 0  # Iteration counter
    # Newton-Raphson iteration to solve for z
    while abs(ratio) > tol and n <= nmax:
        n += 1
        y_val = y(z, r1, r2, A)  # Compute y(z)
        C_val = C(z)             # Compute Stumpff C(z)
        # Check for valid values to avoid sqrt of negative numbers
        if y_val <= 0 or C_val <= 0:
            print(f"Warning: y({z})={y_val}, C({z})={C_val} not positive. Breaking iteration.")
            break
        dFdz_val = dFdz(z, r1, r2, A, mu)  # Derivative of F(z)
        F_val = F(z, t, r1, r2, A, mu)     # Value of F(z)
        # If values are invalid, break
        if np.isnan(F_val) or np.isnan(dFdz_val) or dFdz_val == 0:
            print(f"Warning: F({z})={F_val}, dFdz({z})={dFdz_val} invalid. Breaking iteration.")
            break
        ratio = F_val / dFdz_val  # Newton's update
        z -= ratio                # Update z

    # Warn if Newton's method did not converge
    if n >= nmax:
        print('\n\n **Number of iterations exceeds {} in lambert \n\n '.format(nmax))

    # Compute Lagrange coefficients (f, g, gdot)
    # These relate position and velocity vectors for the transfer
    f = 1 - y(z, r1, r2, A) / r1
    g = A * np.sqrt(y(z, r1, r2, A) / mu)
    gdot = 1 - y(z, r1, r2, A) / r2

    # Compute initial and final velocity vectors using Lagrange coefficients
    V1 = (R2 - f * R1) / g  # Initial velocity vector
    V2 = (gdot * R2 - R1) / g  # Final velocity vector

    # Compute true anomalies for R1 and R2
    def true_anomaly(R):
        # Returns true anomaly in radians for position vector R
        return np.arctan2(R[1], R[0])
    theta1 = true_anomaly(R1)
    theta2 = true_anomaly(R2)

    return V1, V2, theta1, theta2, R1, R2  # Return velocity vectors and true anomalies

def y(z, r1, r2, A):
    """
    Computes y(z) as a function of universal anomaly z.
    Equation 5.38 from Bate, Mueller, White.

    Parameters:
        z (float): Universal anomaly
        r1, r2 (float): Magnitudes of position vectors
        A (float): Geometry constant
    Returns:
        y (float): Value used in universal variable formulation
    """
    # y(z) combines geometry and Stumpff functions for the transfer
    return r1 + r2 + A * (z * S(z) - 1) / np.sqrt(C(z))

def F(z, t, r1, r2, A, mu):
    """
    Computes F(z, t) for Newton's method.
    Equation 5.40 from Bate, Mueller, White.

    Parameters:
        z (float): Universal anomaly
        t (float): Time of flight
        r1, r2, A, mu: As above
    Returns:
        F (float): Value for Newton's method root finding
    """
    # F(z, t) is the function whose root gives the correct z for the transfer
    return (y(z, r1, r2, A) / C(z)) ** 1.5 * S(z) + A * np.sqrt(y(z, r1, r2, A)) - np.sqrt(mu) * t

def dFdz(z, r1, r2, A, mu):
    """
    Computes derivative dF/dz for Newton's method.
    Equation 5.43 from Bate, Mueller, White.

    Parameters:
        z (float): Universal anomaly
        r1, r2, A, mu: As above
    Returns:
        dFdz (float): Derivative for Newton's method
    """
    # Special case for z = 0 (series expansion)
    if z == 0:
        return np.sqrt(2) / 40 * y(0, r1, r2, A) ** 1.5 + \
               A / 8 * (np.sqrt(y(0, r1, r2, A)) + A * np.sqrt(1 / 2 / y(0, r1, r2, A)))
    else:
        # General case for z != 0
        return (y(z, r1, r2, A) / C(z)) ** 1.5 * \
               (1 / 2 / z * (C(z) - 3 * S(z) / 2 / C(z)) + 3 * S(z) ** 2 / 4 / C(z)) + \
               A / 8 * (3 * S(z) / C(z) * np.sqrt(y(z, r1, r2, A)) + A * np.sqrt(C(z) / y(z, r1, r2, A)))

def C(z):
    """
    Returns Stumpff function C(z).

    Parameters:
        z (float): Universal anomaly
    Returns:
        C (float): Stumpff function value
    """
    return stumpC(z)

def S(z):
    """
    Returns Stumpff function S(z).

    Parameters:
        z (float): Universal anomaly
    Returns:
        S (float): Stumpff function value
    """
    return stumpS(z)

# --- Stumpff functions ---
def stumpC(z):
    """
    Computes Stumpff function C(z) robustly for all z.
    Handles z ~ 0 with series expansion.

    Parameters:
        z (float): Universal anomaly
    Returns:
        C (float): Stumpff function value
    """
    epsZ = 1e-8  # Small value for series expansion
    if abs(z) < epsZ:
        return 0.5  # Series expansion about z=0
    if z > 0:
        s = np.sqrt(z)
        return (1 - np.cos(s)) / z  # Elliptic case
    else:
        s = np.sqrt(-z)
        return (1 - np.cosh(s)) / z  # Hyperbolic case

def stumpS(z):
    """
    Computes Stumpff function S(z) robustly for all z.
    Handles z ~ 0 with series expansion.

    Parameters:
        z (float): Universal anomaly
    Returns:
        S (float): Stumpff function value
    """
    epsZ = 1e-8  # Small value for series expansion
    if abs(z) < epsZ:
        return 1.0 / 6.0  # Series expansion about z=0
    if z > 0:
        s = np.sqrt(z)
        return (s - np.sin(s)) / (s ** 3)  # Elliptic case
    else:
        s = np.sqrt(-z)
        return (np.sinh(s) - s) / (s ** 3)  # Hyperbolic case


