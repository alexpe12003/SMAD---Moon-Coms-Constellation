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
	tof_days = 10  # float(input("Enter time of flight in days: "))
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
	R_moon_geo = R_moon_helio - R_earth_helio  # Moon's position relative to Earth's center
	
	# Position parking orbit at the side of Earth relative to Moon direction
	moon_direction = R_moon_geo / np.linalg.norm(R_moon_geo)  # Unit vector towards Moon
	# Create a perpendicular vector in the xy-plane
	perpendicular_direction = np.array([-moon_direction[1], moon_direction[0], 0])
	# Normalize the perpendicular vector
	perpendicular_direction = perpendicular_direction / np.linalg.norm(perpendicular_direction)
	R_earth_geo = R_parking * perpendicular_direction  # Parking orbit at the side of Earth
	
	# Add previous satellite position behind Earth (diametrically opposite to parking orbit)
	R_satellite_previous = -R_parking * moon_direction  # Position behind Earth, opposite to parking orbit

	# Time of flight in seconds
	t_sec = tof_days * 24 * 3600

	# Earth's gravitational parameter (km^3/s^2)
	mu_earth = 3.986004418e5

	print("Geocentric Earth position at departure (km):", R_earth_geo)
	print("Geocentric Moon position at arrival (km):", R_moon_geo)
	print("Previous satellite position behind Earth (km):", R_satellite_previous)

	V1_geo, V2_geo, theta1_lambert, theta2_lambert = lambert(R_earth_geo, R_moon_geo, t_sec, string='pro', mu=mu_earth)
	print("\nGeocentric Lambert solution (parking orbit):")
	print("Initial velocity (km/s):", V1_geo)
	print("Initial velocity magnitude (km/s):", np.linalg.norm(V1_geo))
	print("Final velocity (km/s):", V2_geo)
	print("Final velocity magnitude (km/s):", np.linalg.norm(V2_geo))
	
	# Calculate delta-V requirements
	# 1. Circular parking orbit velocity at departure
	V_parking_circular = np.sqrt(mu_earth / R_parking)  # Circular velocity at parking orbit
	V_parking_vector = V_parking_circular * np.array([-perpendicular_direction[1], perpendicular_direction[0], 0])  # Tangential velocity
	
	# 2. Delta-V at departure (parking orbit to transfer orbit)
	delta_V1 = V1_geo - V_parking_vector
	delta_V1_magnitude = np.linalg.norm(delta_V1)
	
	# 3. Moon's orbital velocity (for reference at arrival)
	# Assuming Moon is in circular orbit around Earth (approximation)
	R_moon_distance = np.linalg.norm(R_moon_geo)
	V_moon_circular = np.sqrt(mu_earth / R_moon_distance)
	# Moon velocity vector (perpendicular to position vector)
	moon_unit = R_moon_geo / R_moon_distance
	V_moon_vector = V_moon_circular * np.array([-moon_unit[1], moon_unit[0], 0])
	
	# 4. Delta-V at arrival (transfer orbit to Moon orbit)
	delta_V2 = V2_geo - V_moon_vector
	delta_V2_magnitude = np.linalg.norm(delta_V2)
	
	# Total delta-V
	total_delta_V = delta_V1_magnitude + delta_V2_magnitude
	
	print("\n=== DELTA-V ANALYSIS ===")
	print(f"Parking orbit circular velocity: {V_parking_circular:.3f} km/s")
	print(f"Parking orbit velocity vector: [{V_parking_vector[0]:.3f}, {V_parking_vector[1]:.3f}, {V_parking_vector[2]:.3f}] km/s")
	print(f"Moon orbital velocity (circular approx): {V_moon_circular:.3f} km/s")
	print(f"Moon velocity vector: [{V_moon_vector[0]:.3f}, {V_moon_vector[1]:.3f}, {V_moon_vector[2]:.3f}] km/s")
	print(f"\nDelta-V at departure: {delta_V1_magnitude:.3f} km/s")
	print(f"Delta-V vector at departure: [{delta_V1[0]:.3f}, {delta_V1[1]:.3f}, {delta_V1[2]:.3f}] km/s")
	print(f"Delta-V at arrival: {delta_V2_magnitude:.3f} km/s")
	print(f"Delta-V vector at arrival: [{delta_V2[0]:.3f}, {delta_V2[1]:.3f}, {delta_V2[2]:.3f}] km/s")
	print(f"\nTOTAL DELTA-V: {total_delta_V:.3f} km/s")

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

	# Plot Geocentric Earth and Moon positions
	fig = plt.figure(figsize=(10, 8))
	ax = fig.add_subplot(111, projection='3d')
	
	# Plot Earth at origin
	ax.scatter(0, 0, 0, color='blue', s=100, label='Earth (Center)')
	
	# Plot previous satellite position behind Earth
	ax.scatter(R_satellite_previous[0], R_satellite_previous[1], R_satellite_previous[2], 
	          color='red', s=40, label='Satellite (Previous Position)')
	
	# Plot parking orbit position
	ax.scatter(R_earth_geo[0], R_earth_geo[1], R_earth_geo[2], color='green', s=50, label='Parking Orbit (Start)')
	
	# Plot Moon position
	ax.scatter(R_moon_geo[0], R_moon_geo[1], R_moon_geo[2], color='gray', s=80, label='Moon (Arrival)')
	
	# Plot vectors from Earth center
	ax.plot([0, R_satellite_previous[0]], [0, R_satellite_previous[1]], [0, R_satellite_previous[2]], 
	        'r--', alpha=0.7, label='Earth to Previous Position')
	ax.plot([0, R_earth_geo[0]], [0, R_earth_geo[1]], [0, R_earth_geo[2]], 'g--', alpha=0.7, label='Earth to Parking Orbit')
	ax.plot([0, R_moon_geo[0]], [0, R_moon_geo[1]], [0, R_moon_geo[2]], 'gray', alpha=0.7, label='Earth to Moon')
	
	# Plot satellite trajectory from previous position to parking orbit
	ax.plot([R_satellite_previous[0], R_earth_geo[0]], [R_satellite_previous[1], R_earth_geo[1]], 
	        [R_satellite_previous[2], R_earth_geo[2]], 'orange', linewidth=2, alpha=0.8, label='Satellite Path')
	
	ax.set_xlabel('X (km)')
	ax.set_ylabel('Y (km)')
	ax.set_zlabel('Z (km)')
	ax.set_title('Geocentric Earth-Moon Transfer')
	ax.legend()
	plt.tight_layout()
	plt.show()


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
