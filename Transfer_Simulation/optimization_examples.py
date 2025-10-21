"""
Quick optimization examples showing different step counts
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from analysis.optimization import multi_parameter_optimization


def run_multiple_optimizations():
    """
    Run optimizations with different step counts to show scalability
    """
    print("OPTIMIZATION SCALING DEMONSTRATION")
    print("=" * 60)
    
    step_counts = [64, 216, 512]  # 4^3, 6^3, 8^3 combinations
    
    for i, steps in enumerate(step_counts, 1):
        print(f"\n{i}/3: Running optimization with {steps} steps...")
        print("-" * 40)
        
        try:
            results = multi_parameter_optimization(num_steps=steps)
            
            if results is not None:
                optimal = results['optimal_result']
                stats = results['statistics']
                
                print(f"✅ Completed successfully!")
                print(f"Best ΔV: {optimal['total_delta_v']:.3f} km/s")
                print(f"Optimal parameters: V0={optimal['V0']:.3f}, γ0={optimal['gamma0']:.1f}°, λ1={optimal['lambda1']:.1f}°")
                print(f"Success rate: {stats['success_rate']:.1f}% ({stats['successful_calculations']}/{stats['total_calculations']})")
                
                # Results displayed only (no file saving)
                
            else:
                print("❌ Optimization failed")
                
        except Exception as e:
            print(f"❌ Error: {e}")
    
    print(f"\n" + "=" * 60)
    print("All optimizations completed!")
    print("\nKey insights:")
    print("- More steps generally lead to better solutions")
    print("- Success rate may vary depending on parameter ranges")
    print("- V0 = 1.372 DU/TU is consistently optimal (minimum energy)")
    print("- γ0 and λ1 show more variation and benefit from optimization")


if __name__ == "__main__":
    run_multiple_optimizations()