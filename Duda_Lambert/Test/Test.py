import numpy as np

# Orbit A parameters
e = 0.1  # eccentricity
a = 7000.0  # km, semi-major axis
mu = 398600.4418  # km^3/s^2, Earth's gravitational parameter

# True anomalies (degrees)
theta1_deg = 30.0
theta2_deg = 60.0

theta1 = np.radians(theta1_deg)
theta2 = np.radians(theta2_deg)

# a) Radii at theta1 and theta2
r1 = a * (1 - e**2) / (1 + e * np.cos(theta1))
r2 = a * (1 - e**2) / (1 + e * np.cos(theta2))

# b) Velocities at theta1 and theta2 (vis-viva equation)
v1 = np.sqrt(mu * (2/r1 - 1/a))
v2 = np.sqrt(mu * (2/r2 - 1/a))

# c) Time of flight between theta1 and theta2
# First, compute eccentric anomaly E1 and E2
cos_E1 = (e + np.cos(theta1)) / (1 + e * np.cos(theta1))
E1 = np.arccos(cos_E1)
if theta1 > np.pi:
    E1 = 2 * np.pi - E1

cos_E2 = (e + np.cos(theta2)) / (1 + e * np.cos(theta2))
E2 = np.arccos(cos_E2)
if theta2 > np.pi:
    E2 = 2 * np.pi - E2

# Mean anomalies
M1 = E1 - e * np.sin(E1)
M2 = E2 - e * np.sin(E2)

# Time of flight (in seconds)
T = 2 * np.pi * np.sqrt(a**3 / mu)
tof = (M2 - M1) * T / (2 * np.pi)
if tof < 0:
    tof += T

# Print results

print(f"Orbit A: e = {e}, a = {a} km")
print(f"True anomalies: theta1 = {theta1_deg} deg, theta2 = {theta2_deg} deg")
print(f"Initial radius r1 = {r1:.2f} km")
print(f"Final radius r2 = {r2:.2f} km")
print(f"Initial velocity v1 = {v1:.2f} km/s")
print(f"Final velocity v2 = {v2:.2f} km/s")
print(f"Time of flight (tof) = {tof:.2f} s")

# Lambert check
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../lambert')))
from lambert import lambert


# Construct position vectors (in plane, arbitrary orientation)
R1 = np.array([r1 * np.cos(theta1), r1 * np.sin(theta1), 0])
R2 = np.array([r2 * np.cos(theta2), r2 * np.sin(theta2), 0])

V1_lambert, V2_lambert, theta1_lambert, theta2_lambert = lambert(R1, R2, tof, string='pro', mu=mu)

print("\nLambert solution:")
print(f"Initial velocity (km/s): {np.linalg.norm(V1_lambert):.2f}")
print(f"Final velocity (km/s): {np.linalg.norm(V2_lambert):.2f}")
print(f"Initial theta (deg): {np.degrees(theta1_lambert):.2f}")
print(f"Final theta (deg): {np.degrees(theta2_lambert):.2f}")

# Calculate classical orbital elements from Lambert initial state
aux_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../Auxiliary'))
if aux_path not in sys.path:
    sys.path.append(aux_path)
from coe_sv import coe_from_sv

coe_lambert = coe_from_sv(R1, V1_lambert, mu)
print(f"Lambert Orbit: e: {coe_lambert[1]:.4f}, a (km): {coe_lambert[6]:.2f}, ")