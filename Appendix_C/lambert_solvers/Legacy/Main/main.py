"""
Main script to solve the Lambert problem for an Earth-Moon transfer.
Uses Skyfield to get ephemeris data for Earth and Moon.
"""


from skyfield.api import load, Topos
from datetime import datetime, timedelta
import numpy as np

import sys
import os
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

if __name__ == "__main__":
	main()
