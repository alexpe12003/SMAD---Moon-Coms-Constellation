"""
Quick RK4 vs RK8 Comparison Plot
================================

Creates a simple but powerful comparison plot showing the superior
accuracy of RK8 over RK4 for orbital mechanics integration.
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.Constants import INTEGRATION_STEPS, RESULTS_DIR
from main_simulation import OrbitSimulation

# =============================================================================
# QUICK COMPARISON PLOT
# =============================================================================

def create_rk4_rk8_comparison():
    """Create a focused comparison plot of RK4 vs RK8 accuracy and performance."""
    
    # Load results from a recent simulation
    sim = OrbitSimulation()
    sim.setup_initial_conditions()
    convergence_results = sim.run_convergence_study()
    
    step_sizes = np.array(INTEGRATION_STEPS)
    
    # Extract data for both methods
    rk4_pos_errors = [convergence_results['RK4'][h]['max_position_error'] for h in step_sizes]
    rk8_pos_errors = [convergence_results['RK8'][h]['max_position_error'] for h in step_sizes]
    rk4_times = [convergence_results['RK4'][h]['computation_time'] for h in step_sizes]
    rk8_times = [convergence_results['RK8'][h]['computation_time'] for h in step_sizes]
    
    # Create comparison plots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # 1. Accuracy Comparison
    ax1.loglog(step_sizes, rk4_pos_errors, 'ro-', label='RK4', linewidth=2, markersize=8)
    ax1.loglog(step_sizes, rk8_pos_errors, 'bo-', label='RK8', linewidth=2, markersize=8)
    ax1.set_xlabel('Step Size (s)', fontsize=12)
    ax1.set_ylabel('Maximum Position Error (km)', fontsize=12)
    ax1.set_title('Accuracy Comparison: RK4 vs RK8', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=12)
    ax1.grid(True, alpha=0.3)
    
    # 2. Computation Time vs Step Size
    ax2.loglog(step_sizes, rk4_times, 'ro-', label='RK4', linewidth=2, markersize=8)
    ax2.loglog(step_sizes, rk8_times, 'bo-', label='RK8', linewidth=2, markersize=8)
    ax2.set_xlabel('Step Size (s)', fontsize=12)
    ax2.set_ylabel('Computation Time (s)', fontsize=12)
    ax2.set_title('Computation Time vs Step Size', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=12)
    ax2.grid(True, alpha=0.3)
    
    # 3. Efficiency Plot: Accuracy vs Time
    ax3.loglog(rk4_times, rk4_pos_errors, 'ro-', label='RK4', linewidth=2, markersize=8)
    ax3.loglog(rk8_times, rk8_pos_errors, 'bo-', label='RK8', linewidth=2, markersize=8)
    ax3.set_xlabel('Computation Time (s)', fontsize=12)
    ax3.set_ylabel('Maximum Position Error (km)', fontsize=12)
    ax3.set_title('Efficiency: Accuracy vs Computation Time', fontsize=14, fontweight='bold')
    ax3.legend(fontsize=12)
    ax3.grid(True, alpha=0.3)
    
    # 4. Relative Performance Table
    ax4.axis('off')
    
    # Create performance summary
    table_data = []
    table_data.append(['Step Size (s)', 'RK4 Error (km)', 'RK8 Error (km)', 'RK8/RK4 Accuracy Ratio'])
    
    for i, h in enumerate([1.0, 10.0, 100.0, 300.0]):
        if h in step_sizes:
            rk4_err = convergence_results['RK4'][h]['max_position_error']
            rk8_err = convergence_results['RK8'][h]['max_position_error']
            ratio = rk4_err / rk8_err
            table_data.append([f'{h:g}', f'{rk4_err:.2e}', f'{rk8_err:.2e}', f'{ratio:.1f}x'])
    
    # Create table
    table = ax4.table(cellText=table_data[1:], colLabels=table_data[0],
                     cellLoc='center', loc='center',
                     colWidths=[0.2, 0.25, 0.25, 0.3])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)
    
    # Style the header
    for i in range(len(table_data[0])):
        table[(0, i)].set_facecolor('#40466e')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    ax4.set_title('Performance Summary', fontsize=14, fontweight='bold', pad=20)
    
    plt.tight_layout()
    plt.savefig(f'{RESULTS_DIR}/rk4_vs_rk8_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Print summary
    print("\n" + "="*80)
    print("RK4 vs RK8 COMPARISON SUMMARY")
    print("="*80)
    print(f"For the 300s step size:")
    rk4_300_err = convergence_results['RK4'][300.0]['max_position_error']
    rk8_300_err = convergence_results['RK8'][300.0]['max_position_error']
    rk4_300_time = convergence_results['RK4'][300.0]['computation_time']
    rk8_300_time = convergence_results['RK8'][300.0]['computation_time']
    
    print(f"  RK4: {rk4_300_err:.1e} km error in {rk4_300_time:.3f} s")
    print(f"  RK8: {rk8_300_err:.1e} km error in {rk8_300_time:.3f} s")
    print(f"  RK8 is {rk4_300_err/rk8_300_err:.0f}x more accurate")
    print(f"  RK8 is {rk8_300_time/rk4_300_time:.1f}x slower")
    print(f"\nConclusion: RK8 provides dramatically better accuracy for larger step sizes")
    print("="*80)

if __name__ == "__main__":
    # Ensure results directory exists
    os.makedirs(RESULTS_DIR, exist_ok=True)
    
    print("Creating RK4 vs RK8 comparison...")
    create_rk4_rk8_comparison()
    print(f"Comparison plot saved to {RESULTS_DIR}/rk4_vs_rk8_comparison.png")