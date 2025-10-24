#!/usr/bin/env python3
"""
Test the modified optimization function with fixed V0 value
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analysis.optimization import multi_parameter_optimization

def test_fixed_v0_optimization():
    """
    Test the optimization with V0 fixed at 1.372 DU/TU
    """
    print("Testing optimization with V0 fixed at 1.372 DU/TU")
    print("=" * 60)
    
    # Run optimization with a smaller number of steps for testing
    results = multi_parameter_optimization(num_steps=100)  # This will be sqrt(100) = 10 steps per parameter
    
    if results is None:
        print("ERROR: Optimization failed!")
        return
    
    print("\nTest completed successfully!")
    
    optimal = results['optimal_result']
    print(f"\nOptimal parameters found:")
    print(f"  R0 = {optimal['R0']:.3f} DU (fixed)")
    print(f"  V0 = {optimal['V0']:.3f} DU/TU (fixed)")
    print(f"  γ0 = {optimal['gamma0']:.3f}°")
    print(f"  λ1 = {optimal['lambda1']:.3f}°")
    print(f"  Total ΔV = {optimal['total_delta_v']:.3f} km/s")
    
    # Verify V0 is actually fixed
    all_v0_values = [r['V0'] for r in results['all_results']]
    unique_v0_values = set(all_v0_values)
    
    print(f"\nVerification:")
    print(f"  All V0 values are identical: {len(unique_v0_values) == 1}")
    print(f"  V0 value used: {list(unique_v0_values)}")
    print(f"  Total successful calculations: {results['statistics']['successful_calculations']}")

if __name__ == "__main__":
    test_fixed_v0_optimization()