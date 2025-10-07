"""
Configuration file containing constants and mission parameters for lunar transfer analysis
"""

import math

# Physical constants
EARTH_RADIUS_KM = 6378  # km
MOON_RADIUS_KM = 1737   # km
MU_EARTH = 1.0          # Earth gravitational parameter in canonical units
MU_MOON_KMS = 4902.829028  # Moon gravitational parameter in km³/s²

# Lunar orbit constants
MOON_DISTANCE = 60.27   # DU - Distance from Earth to Moon
MOON_SOI_RADIUS = 10.395  # DU - Moon's sphere of influence radius
MOON_ORBITAL_VELOCITY = 1.018  # km/s - Moon's orbital velocity relative to Earth

# Mission parameters
DEFAULT_PARKING_ORBIT_RADIUS = 1.05  # DU
DEFAULT_TRANSFER_VELOCITY = 1.372    # DU/TU
DEFAULT_FLIGHT_PATH_ANGLE = 0        # degrees
DEFAULT_TARGET_PERIGEE_ALTITUDE = 600  # km

# Conversion factors
DU_TO_KM = EARTH_RADIUS_KM  # 1 DU = 6378 km
TU_TO_SEC = 806.81          # 1 TU = 806.81 seconds
TU_TO_HOURS = TU_TO_SEC / 3600  # 1 TU in hours
DU_TU_TO_KM_S = DU_TO_KM / TU_TO_SEC  # DU/TU to km/s conversion

# Analysis parameters
LAMBDA1_RANGE = (0, 360)    # degrees
LAMBDA1_STEP = 5            # degrees