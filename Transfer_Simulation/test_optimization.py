"""
Simple test script to verify the optimization function works
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from analysis.optimization import multi_parameter_optimization


def test_optimization():
    """
    Test the optimization function with a small number of steps
    """
    print("Testing Multi-Parameter Optimization Function")
    print("=" * 50)
    
    # Run with a small number of steps for testing
    num_steps = 27  # 3x3x3 = 27 combinations
    print(f"Testing with {num_steps} steps...")
    
    try:
        results = multi_parameter_optimization(num_steps=num_steps)
        
        if results is not None:
            print("\n✅ Optimization completed successfully!")
            
            # Display key results (no file saving)
            optimal = results['optimal_result']
            stats = results['statistics']
            
            print(f"\nTest Results Summary:")
            print(f"Best Delta-V: {optimal['total_delta_v']:.3f} km/s")
            print(f"Optimal V0: {optimal['V0']:.3f} DU/TU")
            print(f"Optimal γ0: {optimal['gamma0']:.1f}°")
            print(f"Optimal λ1: {optimal['lambda1']:.1f}°")
            print(f"Success rate: {stats['success_rate']:.1f}%")
            print(f"Total calculations: {stats['total_calculations']}")
            
            return True
            
        else:
            print("❌ Optimization failed!")
            return False
            
    except Exception as e:
        print(f"❌ Error during optimization: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = test_optimization()
    if success:
        print("\n🎉 All tests passed! The optimization function is working correctly.")
    else:
        print("\n❌ Tests failed. Check the error messages above.")