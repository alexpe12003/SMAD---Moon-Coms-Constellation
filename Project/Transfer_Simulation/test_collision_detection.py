#!/usr/bin/env python3
"""
Test the improved collision detection in the optimization function
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analysis.optimization import multi_parameter_optimization

def test_improved_collision_detection():
    """
    Test the optimization with improved collision detection
    """
    print("Testing improved collision detection in optimization...")
    print("=" * 60)
    
    # Run optimization with small number of steps for testing
    results = multi_parameter_optimization(num_steps=400)  # 20x20 grid
    
    if results is None:
        print("ERROR: Optimization failed!")
        return
    
    print("\n" + "=" * 60)
    print("COLLISION DETECTION TEST RESULTS")
    print("=" * 60)
    
    stats = results['statistics']
    optimal = results['optimal_result']
    
    print(f"Total calculations: {stats['total_calculations']}")
    print(f"Successful: {stats['successful_calculations']} ({stats['success_rate']:.1f}%)")
    print(f"Collision trajectories: {stats['collision_count']} ({stats['collision_rate']:.1f}%)")
    print(f"Invalid trajectories: {stats['invalid_count']} ({stats['invalid_rate']:.1f}%)")
    print()
    
    print("OPTIMAL RESULT VALIDATION:")
    print(f"  Parameters: γ0={optimal['gamma0']:.1f}°, λ1={optimal['lambda1']:.1f}°")
    print(f"  Total ΔV: {optimal['total_delta_v']:.3f} km/s")
    print(f"  Hyperbolic perigee: {optimal.get('rp_hyperbolic_km', 'N/A'):.1f} km")
    print(f"  Safety margin: {optimal.get('rp_hyperbolic_km', 0) - 1737.0:.1f} km")
    print(f"  Orbit type: {optimal.get('orbit_type', 'unknown')}")
    
    # Validation checks
    rp_hyp = optimal.get('rp_hyperbolic_km', 0)
    if rp_hyp > 1737.0:
        print("✅ PASS: Optimal result is collision-free")
    else:
        print("❌ FAIL: Optimal result has collision trajectory!")
        
    if optimal['total_delta_v'] < 10.0:
        print("✅ PASS: Optimal delta-V is reasonable")
    else:
        print("❌ FAIL: Optimal delta-V is too high!")
    
    return results

if __name__ == "__main__":
    test_improved_collision_detection()