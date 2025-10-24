"""
Simple Porkchop Plot Test
========================
A simplified version to test the concept and verify Lambert solver integration.
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import sys
import os

# Add the lambert solver to the path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'Duda_Lambert', 'lambert')))

# Test imports
try:
    import lambert
    from porkchop import create_leo_satellite, create_moon_orbit
    print("All imports successful!")
except ImportError as e:
    print(f"Import error: {e}")
    sys.exit(1)

# Physical Constants
MU_EARTH_KM3_S2 = 398600.4418  # Earth's gravitational parameter (km³/s²)

def test_lambert_solver():
    """Test the Lambert solver with a simple example"""
    print("\nTesting Lambert solver...")
    
    # Simple test case - transfer from LEO to a higher orbit
    R1 = np.array([7000, 0, 0])      # Initial position (km)
    R2 = np.array([0, 10000, 0])     # Final position (km)
    tof = 3600                        # Time of flight (s) - 1 hour
    
    try:
        V1, V2, theta1, theta2 = lambert.lambert(R1, R2, tof, 'pro', MU_EARTH_KM3_S2)
        print(f"Lambert solver test successful!")
        print(f"V1 magnitude: {np.linalg.norm(V1):.3f} km/s")
        print(f"V2 magnitude: {np.linalg.norm(V2):.3f} km/s")
        return True
    except Exception as e:
        print(f"Lambert solver test failed: {e}")
        return False

def simple_porkchop_test():
    """Generate a simple porkchop plot with reduced complexity"""
    print("\nGenerating simple porkchop test...")
    
    # Create orbit objects
    leo_satellite = create_leo_satellite()
    moon = create_moon_orbit()
    
    # Simple parameters
    departure_times = np.linspace(0, leo_satellite.period, 10)  # 10 points over one LEO orbit
    tof_range = np.linspace(24*3600, 120*3600, 10)            # 1 to 5 days
    
    print(f"Departure times: {len(departure_times)} points")
    print(f"Time of flight: {len(tof_range)} points")
    
    delta_v_matrix = np.zeros((len(tof_range), len(departure_times)))
    
    for i, dep_time in enumerate(departure_times):
        for j, tof in enumerate(tof_range):
            # Get LEO position at departure
            theta_rad, _, x1, y1 = leo_satellite.get_position_at_time(dep_time)
            R1 = np.array([x1, y1, 0])
            
            # Get Moon position at arrival
            arrival_time = dep_time + tof
            theta_rad, _, x2, y2 = moon.get_position_at_time(arrival_time)
            R2 = np.array([x2, y2, 0])
            
            try:
                # Solve Lambert's problem
                V1, V2, _, _ = lambert.lambert(R1, R2, tof, 'pro', MU_EARTH_KM3_S2)
                
                # Calculate delta-V
                r1 = np.linalg.norm(R1)
                r2 = np.linalg.norm(R2)
                
                v1_circular = math.sqrt(MU_EARTH_KM3_S2 / r1)
                v2_circular = math.sqrt(MU_EARTH_KM3_S2 / r2)
                
                v1_transfer = np.linalg.norm(V1)
                v2_transfer = np.linalg.norm(V2)
                
                delta_v1 = abs(v1_transfer - v1_circular)
                delta_v2 = abs(v2_transfer - v2_circular)
                delta_v_total = delta_v1 + delta_v2
                
                delta_v_matrix[j, i] = delta_v_total
                
            except Exception as e:
                print(f"Error at dep_time={dep_time/60:.1f}min, tof={tof/3600:.1f}h: {e}")
                delta_v_matrix[j, i] = np.nan
    
    # Create plot
    plt.figure(figsize=(10, 6))
    
    departure_grid, tof_grid = np.meshgrid(departure_times/60, tof_range/3600)  # Convert to minutes and hours
    
    # Filter out invalid values
    valid_mask = ~np.isnan(delta_v_matrix) & (delta_v_matrix < 20)  # Reasonable delta-V limit
    
    if np.any(valid_mask):
        # Create contour plot
        contour = plt.contourf(departure_grid, tof_grid, delta_v_matrix, 
                             levels=20, cmap='viridis')
        plt.colorbar(contour, label='Total ΔV (km/s)')
        
        # Add contour lines
        plt.contour(departure_grid, tof_grid, delta_v_matrix, 
                   levels=20, colors='white', alpha=0.5, linewidths=0.5)
        
        # Find minimum
        min_idx = np.unravel_index(np.nanargmin(delta_v_matrix), delta_v_matrix.shape)
        min_dep = departure_grid[min_idx]
        min_tof = tof_grid[min_idx]
        min_dv = delta_v_matrix[min_idx]
        
        plt.plot(min_dep, min_tof, 'r*', markersize=15, label=f'Min ΔV: {min_dv:.3f} km/s')
        
        print(f"\nOptimal solution found:")
        print(f"Departure: {min_dep:.1f} minutes")
        print(f"Time of flight: {min_tof:.1f} hours")
        print(f"Total ΔV: {min_dv:.3f} km/s")
        
    else:
        print("No valid solutions found")
    
    plt.xlabel('Departure Time (minutes)')
    plt.ylabel('Time of Flight (hours)')
    plt.title('Simple LEO-to-Moon Porkchop Plot')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()
    
    return delta_v_matrix

def main():
    """Main test function"""
    print("SIMPLE PORKCHOP PLOT TEST")
    print("=" * 30)
    
    # Test Lambert solver first
    if not test_lambert_solver():
        print("Lambert solver test failed - cannot continue")
        return
    
    # Generate simple porkchop plot
    result = simple_porkchop_test()
    
    print("\nTest complete!")
    return result

if __name__ == "__main__":
    result = main()