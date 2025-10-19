"""
Orbital Mechanics Core Module
============================

This module contains the core numerical integration and orbital mechanics functions
for the SMAD Moon Communications Constellation study.

Core Components:
- Constants: Physical constants and orbital parameters
- RK4: 4th-order Runge-Kutta numerical integration
- RK8: 8th-order Runge-Kutta (Dormand-Prince) integration  
- Kepler: Analytical orbital mechanics solutions
- coe_sv: Robust orbital elements conversion (Algorithm 4.1)
- analysis_utils: Analysis and plotting utilities

Usage:
    from core.Constants import MU_EARTH, get_orbit_b_parameters
    from core.RK4 import propagate_orbit_rk4
    from core.RK8 import propagate_orbit_rk8
    from core.Kepler import analytical_propagation
"""

# Import core modules for easy access
from .Constants import *
from .RK4 import *
from .RK8 import *
from .Kepler import *
from .coe_sv import *
from .analysis_utils import *

__version__ = "1.0.0"
__author__ = "SMAD Moon Communications Constellation Study"
__description__ = "High-precision orbital mechanics and numerical integration"