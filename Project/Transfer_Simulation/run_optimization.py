"""
Quick optimization script to run the multi-parameter optimization directly
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from analysis.optimization import multi_parameter_optimization


def run_optimization():
    """
    Run the multi-parameter optimization with default settings
    """
    print("Running Multi-Parameter Optimization")
    print("=" * 50)
    print("This will optimize V0, gamma0, and lambda1 to minimize mission delta-V")
    print("R0 is fixed at 1.05 DU")
    print()
    
    # Get number of steps from user or use default
    try:
        num_steps = int(input("Enter number of optimization steps (default 500): ") or "500")
    except ValueError:
        num_steps = 500
        print("Using default value: 500 steps")
    
    # Run optimization
    print(f"\nRunning optimization with {num_steps} steps...")
    results = multi_parameter_optimization(num_steps=num_steps)
    
    # Display results (no file saving)
    if results is not None:
        print("\nOptimization completed successfully!")
        
        # Display key results
        optimal = results['optimal_result']
        print(f"\nBest Delta-V found: {optimal['total_delta_v']:.3f} km/s")
        print(f"Optimal parameters:")
        print(f"  V0 = {optimal['V0']:.3f} DU/TU")
        print(f"  γ0 = {optimal['gamma0']:.1f}°")
        print(f"  λ1 = {optimal['lambda1']:.1f}°")
        
    else:
        print("Optimization failed!")
    
    return results


if __name__ == "__main__":
    results = run_optimization()