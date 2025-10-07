"""
Plotting and visualization functions for lunar transfer analysis
"""

import matplotlib.pyplot as plt
import numpy as np


def create_parametric_plots(lambda1_values, delta_v_values, tof_values):
    """
    Create plots showing the relationship between lambda1 and complete mission parameters
    
    Parameters:
    - lambda1_values: List of lambda1 values in degrees
    - delta_v_values: List of complete mission delta-V values in km/s
    - tof_values: List of complete mission time values in hours
    """
    # Create figure with subplots
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 15))
    
    # Plot 1: Lambda1 vs Complete Mission Delta-V
    ax1.plot(lambda1_values, delta_v_values, 'b-', linewidth=2, marker='o', markersize=3)
    ax1.set_xlabel('Lambda1 (degrees)')
    ax1.set_ylabel('Complete Mission Delta-V (km/s)')
    ax1.set_title('Complete Mission Delta-V vs Lambda1')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 360)
    
    # Add minimum delta-V annotation
    min_delta_v_idx = np.argmin(delta_v_values)
    min_delta_v = delta_v_values[min_delta_v_idx]
    min_lambda1 = lambda1_values[min_delta_v_idx]
    ax1.annotate(f'Min ΔV: {min_delta_v:.3f} km/s\nat λ₁ = {min_lambda1}°', 
                xy=(min_lambda1, min_delta_v), xytext=(min_lambda1 + 50, min_delta_v + 0.5),
                arrowprops=dict(arrowstyle='->', color='red'),
                bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.7))
    
    # Plot 2: Lambda1 vs Complete Mission Time
    ax2.plot(lambda1_values, tof_values, 'r-', linewidth=2, marker='s', markersize=3)
    ax2.set_xlabel('Lambda1 (degrees)')
    ax2.set_ylabel('Complete Mission Time (hours)')
    ax2.set_title('Complete Mission Time vs Lambda1')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, 360)
    
    # Add minimum time annotation
    min_tof_idx = np.argmin(tof_values)
    min_tof = tof_values[min_tof_idx]
    min_tof_lambda1 = lambda1_values[min_tof_idx]
    ax2.annotate(f'Min Time: {min_tof:.1f} hrs\nat λ₁ = {min_tof_lambda1}°', 
                xy=(min_tof_lambda1, min_tof), xytext=(min_tof_lambda1 + 50, min_tof + 20),
                arrowprops=dict(arrowstyle='->', color='red'),
                bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.7))
    
    # Plot 3: Delta-V vs Time of Flight (scatter plot)
    scatter = ax3.scatter(tof_values, delta_v_values, c=lambda1_values, cmap='viridis', 
                         s=50, alpha=0.7, edgecolors='black', linewidth=0.5)
    ax3.set_xlabel('Complete Mission Time (hours)')
    ax3.set_ylabel('Complete Mission Delta-V (km/s)')
    ax3.set_title('Mission Trade-off: Delta-V vs Time of Flight')
    ax3.grid(True, alpha=0.3)
    
    # Add colorbar for lambda1 values
    cbar = plt.colorbar(scatter, ax=ax3)
    cbar.set_label('Lambda1 (degrees)')
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig('lambda1_parametric_study.png', dpi=300, bbox_inches='tight')
    print("Plot saved as 'lambda1_parametric_study.png'")
    plt.show()


def create_complete_parametric_plots(lambda1_values, lunar_delta_v_values, complete_delta_v_values, 
                                   earth_moon_time_values, complete_time_values, earth_departure_dv):
    """
    Create comprehensive plots showing complete mission analysis including Earth departure
    
    Parameters:
    - lambda1_values: List of lambda1 values in degrees
    - lunar_delta_v_values: List of lunar-only delta-V values in km/s
    - complete_delta_v_values: List of complete mission delta-V values in km/s
    - earth_moon_time_values: List of Earth-Moon transit times in hours
    - complete_time_values: List of complete mission times in hours
    - earth_departure_dv: Earth departure delta-V in km/s (constant)
    """
    # Create figure with subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # Plot 1: Lambda1 vs Complete Mission Delta-V
    ax1.plot(lambda1_values, complete_delta_v_values, 'b-', linewidth=2, marker='o', markersize=3)
    ax1.axhline(y=earth_departure_dv, color='r', linestyle='--', alpha=0.7, 
                label=f'Earth Departure: {earth_departure_dv:.3f} km/s')
    ax1.set_xlabel('Lambda1 (degrees)')
    ax1.set_ylabel('Complete Mission Delta-V (km/s)')
    ax1.set_title('Complete Mission Delta-V vs Lambda1')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 360)
    ax1.legend()
    
    # Add minimum delta-V annotation
    min_delta_v_idx = np.argmin(complete_delta_v_values)
    min_delta_v = complete_delta_v_values[min_delta_v_idx]
    min_lambda1 = lambda1_values[min_delta_v_idx]
    ax1.annotate(f'Min ΔV: {min_delta_v:.3f} km/s\nat λ₁ = {min_lambda1}°', 
                xy=(min_lambda1, min_delta_v), xytext=(min_lambda1 + 50, min_delta_v + 0.5),
                arrowprops=dict(arrowstyle='->', color='red'),
                bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.7))
    
    # Plot 2: Lambda1 vs Complete Mission Time
    ax2.plot(lambda1_values, complete_time_values, 'g-', linewidth=2, marker='s', markersize=3)
    ax2.set_xlabel('Lambda1 (degrees)')
    ax2.set_ylabel('Complete Mission Time (hours)')
    ax2.set_title('Complete Mission Time vs Lambda1')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, 360)
    
    # Add minimum time annotation
    min_time_idx = np.argmin(complete_time_values)
    min_time = complete_time_values[min_time_idx]
    min_time_lambda1 = lambda1_values[min_time_idx]
    ax2.annotate(f'Min Time: {min_time:.1f} hrs\nat λ₁ = {min_time_lambda1}°', 
                xy=(min_time_lambda1, min_time), xytext=(min_time_lambda1 + 50, min_time + 20),
                arrowprops=dict(arrowstyle='->', color='red'),
                bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.7))
    
    # Plot 3: Lunar-only Delta-V vs Lambda1
    ax3.plot(lambda1_values, lunar_delta_v_values, 'm-', linewidth=2, marker='^', markersize=3)
    ax3.set_xlabel('Lambda1 (degrees)')
    ax3.set_ylabel('Lunar Operations Delta-V (km/s)')
    ax3.set_title('Lunar Operations Delta-V vs Lambda1')
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim(0, 360)
    
    # Plot 4: Delta-V vs Time Trade-off (scatter plot)
    scatter = ax4.scatter(complete_time_values, complete_delta_v_values, c=lambda1_values, 
                         cmap='viridis', s=50, alpha=0.7, edgecolors='black', linewidth=0.5)
    ax4.set_xlabel('Complete Mission Time (hours)')
    ax4.set_ylabel('Complete Mission Delta-V (km/s)')
    ax4.set_title('Mission Trade-off: Delta-V vs Time')
    ax4.grid(True, alpha=0.3)
    
    # Add colorbar for lambda1 values
    cbar = plt.colorbar(scatter, ax=ax4)
    cbar.set_label('Lambda1 (degrees)')
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig('complete_mission_parametric_study.png', dpi=300, bbox_inches='tight')
    print("Complete mission plot saved as 'complete_mission_parametric_study.png'")
    plt.show()