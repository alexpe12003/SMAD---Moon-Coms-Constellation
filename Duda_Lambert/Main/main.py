"""
Main script to solve the Lambert problem for an Earth-Moon transfer.
Uses Skyfield to get ephemeris data for Earth and Moon.
"""


from skyfield.api import load, Topos
from datetime import datetime, timedelta
import numpy as np

import sys
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# Ensure parent directory is in sys.path for package imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../Auxiliary')))

from lambert.lambert import lambert
from Auxiliary.coe_sv import coe_from_sv

def get_body_position(body, ts, eph, date):
	"""
	Returns the position vector of a body (Earth or Moon) in km in the ICRF frame at a given date.
	"""
	t = ts.utc(date.year, date.month, date.day, date.hour, date.minute, date.second)
	# Get position in ICRF (heliocentric)
	pos = eph[body].at(t).position.km
	return np.array(pos)

def main():
	# User input: departure date and time of flight (in days)
	departure_str = "2025-09-21"  # input("Enter departure date (YYYY-MM-DD): ")
	tof_days = 3  # float(input("Enter time of flight in days: "))
	departure_date = datetime.strptime(departure_str, "%Y-%m-%d")
	arrival_date = departure_date + timedelta(days=tof_days)

	# Load ephemeris data
	eph = load('de421.bsp')
	ts = load.timescale()

	# Get heliocentric positions
	R_earth_helio = get_body_position('earth', ts, eph, departure_date)
	R_moon_helio = get_body_position('moon', ts, eph, arrival_date)

	# --- Geocentric transfer with parking orbit ---
	R_earth_radius = 6378.0  # km
	parking_altitude = 200.0  # km (change as needed)
	R_parking = R_earth_radius + parking_altitude
	R_earth_geo = np.array([R_parking, 0, 0.0])  # Parking orbit position (x-axis)
	R_moon_geo = R_moon_helio - R_earth_helio  # Moon's position relative to Earth's center

	# Time of flight in seconds
	t_sec = tof_days * 24 * 3600

	# Earth's gravitational parameter (km^3/s^2)
	mu_earth = 3.986004418e5

	print("Geocentric Earth position at departure (km):", R_earth_geo)
	print("Geocentric Moon position at arrival (km):", R_moon_geo)

	V1_geo, V2_geo = lambert(R_earth_geo, R_moon_geo, t_sec, string='pro', mu=mu_earth)
	print("\nGeocentric Lambert solution (parking orbit):")
	print("Initial velocity (km/s):", V1_geo)
	print("Final velocity (km/s):", V2_geo)

	# Classical Orbital Elements from state vector
	coe = coe_from_sv(R_earth_geo, V1_geo, mu_earth)
	print("\nClassical Orbital Elements (from initial state vector):")
	print("h (km^2/s):", coe[0])
	print("e (eccentricity):", coe[1])
	print("RA (radians):", coe[2])
	print("inclination (radians):", coe[3])
	print("argument of perigee (radians):", coe[4])
	print("true anomaly (radians):", coe[5])
	print("semi-major axis (km):", coe[6])


	# --- Plot geocentric transfer orbit ---
	# Use classical orbital elements to reconstruct the orbit
	# h, e, RA, incl, w, TA, a = coe
	# Generate true anomaly values
	# theta_vals = np.linspace(0, 2 * np.pi, 500)
	# r_vals = a * (1 - e**2) / (1 + e * np.cos(theta_vals))
	# Perifocal coordinates
	# x_pf = r_vals * np.cos(theta_vals)
	# y_pf = r_vals * np.sin(theta_vals)
	# z_pf = np.zeros_like(x_pf)
	# Rotation matrices
	# def rot_matrix(RA, incl, w):
	#     # Rotation from perifocal to geocentric equatorial frame
	#     R3_W = np.array([[np.cos(RA), -np.sin(RA), 0],
	#                     [np.sin(RA),  np.cos(RA), 0],
	#                     [0, 0, 1]])
	#     R1_i = np.array([[1, 0, 0],
	#                     [0, np.cos(incl), -np.sin(incl)],
	#                     [0, np.sin(incl),  np.cos(incl)]])
	#     R3_w = np.array([[np.cos(w), -np.sin(w), 0],
	#                     [np.sin(w),  np.cos(w), 0],
	#                     [0, 0, 1]])
	#     return R3_W @ R1_i @ R3_w
	# Q = rot_matrix(RA, incl, w)
	# Transform to geocentric equatorial frame
	# r_orbit = np.vstack((x_pf, y_pf, z_pf))
	# r_geo = Q @ r_orbit
	# Plot
	# fig = plt.figure(figsize=(8, 6))
	# ax = fig.add_subplot(111, projection='3d')
	# ax.plot(r_geo[0], r_geo[1], r_geo[2], label='Transfer Orbit')
	# ax.scatter(R_earth_geo[0], R_earth_geo[1], R_earth_geo[2], color='blue', label='Parking Orbit (Start)')
	# ax.scatter(R_moon_geo[0], R_moon_geo[1], R_moon_geo[2], color='gray', label='Moon (Arrival)')
	# ax.set_xlabel('X (km)')
	# ax.set_ylabel('Y (km)')
	# ax.set_zlabel('Z (km)')
	# ax.set_title('Geocentric Transfer Orbit (Parking Orbit to Moon)')
	# ax.legend()
	# plt.tight_layout()
	# plt.show()

if __name__ == "__main__":
	main()
