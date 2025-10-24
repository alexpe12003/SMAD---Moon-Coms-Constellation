"""
Test script to verify the new organized structure works correctly
"""

import sys
import os

# Test imports from the new core module
try:
    from core.Constants import MU_EARTH, get_orbit_b_parameters
    print("✅ Constants import successful")
    print(f"   MU_EARTH = {MU_EARTH}")
    
    from core.RK4 import rk4_step, propagate_orbit_rk4
    print("✅ RK4 import successful")
    
    from core.RK8 import rk8_step, propagate_orbit_rk8  
    print("✅ RK8 import successful")
    
    from core.Kepler import analytical_propagation
    print("✅ Kepler import successful")
    
    from core.coe_sv import coe_from_sv
    print("✅ coe_sv import successful")
    
    print("\n🎉 All core module imports successful!")
    print("New folder structure is working correctly!")
    
except ImportError as e:
    print(f"❌ Import error: {e}")
    print("There may be remaining import issues to fix.")