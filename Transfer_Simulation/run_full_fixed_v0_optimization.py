#!/usr/bin/env python3
"""
Run full optimization with V0 fixed at 1.372 DU/TU
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analysis.optimization import multi_parameter_optimization, save_optimization_results

def run_full_optimization():
    """
    Run full optimization with V0 fixed at 1.372 DU/TU and higher resolution
    """
    print("Running FULL optimization with V0 fixed at 1.372 DU/TU")
    print("=" * 70)
    
    # Run optimization with higher resolution 
    # 3600 steps = sqrt(3600) = 60 steps per parameter (gamma0 and lambda1)
    results = multi_parameter_optimization(num_steps=3600)
    
    if results is None:
        print("ERROR: Optimization failed!")
        return
    
    print("\nFull optimization completed successfully!")
    
    # Save results
    save_optimization_results(results, "fixed_v0_optimization_results.npz")
    
    optimal = results['optimal_result']
    print(f"\n" + "=" * 70)
    print("FINAL OPTIMAL PARAMETERS (V0 FIXED)")
    print("=" * 70)
    print(f"  R0 = {optimal['R0']:.3f} DU (fixed)")
    print(f"  V0 = {optimal['V0']:.3f} DU/TU (fixed)")
    print(f"  γ0 = {optimal['gamma0']:.3f}°")
    print(f"  λ1 = {optimal['lambda1']:.3f}°")
    print(f"  Total ΔV = {optimal['total_delta_v']:.3f} km/s")
    print(f"  Earth Departure ΔV = {optimal['earth_departure_dv']:.3f} km/s")
    print(f"  Lunar Operations ΔV = {optimal['lunar_dv']:.3f} km/s")
    print(f"  Total Mission Time = {optimal['total_time_hours']:.1f} hours")

if __name__ == "__main__":
    run_full_optimization()