"""
Visualization tools for the Moon Constellation Deployment System

This script creates plots and diagrams to visualize:
1. Deployment timeline with staggered burns
2. Orbital mechanics of Hohmann transfers  
3. Final constellation geometry
4. Delta-V budget breakdown
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
import os

# Add the satellite deployment module to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from satdep import MoonConstellationDeployer

def plot_deployment_timeline(deployment_schedule, transfer_time_h, stagger_time_h):
    """
    Create a timeline plot showing staggered satellite deployments
    """
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Color palette for satellites
    colors = plt.cm.tab10(np.linspace(0, 1, len(deployment_schedule)))
    
    # Plot deployment timeline for each satellite
    for i, sat in enumerate(deployment_schedule):
        sat_id = sat['satellite_id']
        start_time = sat['deployment_time_h']
        burn1_time = start_time
        burn2_time = start_time + transfer_time_h
        
        # Timeline bar
        ax.barh(i, burn2_time, left=0, height=0.6, color=colors[i], alpha=0.3, 
               label=f'Satellite {sat_id}')
        
        # Burn markers
        ax.scatter(burn1_time, i, color=colors[i], s=100, marker='o', 
                  label='Burn 1' if i == 0 else "")
        ax.scatter(burn2_time, i, color=colors[i], s=100, marker='s', 
                  label='Burn 2' if i == 0 else "")
        
        # Text annotations
        ax.text(burn1_time - 1, i, f'Start', ha='right', va='center', fontsize=8)
        ax.text(burn2_time + 1, i, f'Complete\n({sat["final_separation_deg"]:.0f}°)', 
               ha='left', va='center', fontsize=8)
    
    # Formatting
    ax.set_yticks(range(len(deployment_schedule)))
    ax.set_yticklabels([f'Sat {i+1}' for i in range(len(deployment_schedule))])
    ax.set_xlabel('Time (hours)')
    ax.set_title('Moon Constellation Deployment Timeline\nStaggered Hohmann Transfers', 
                fontsize=14, pad=20)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right')
    
    # Add stagger time annotation
    ax.annotate('', xy=(stagger_time_h, 1), xytext=(0, 1),
               arrowprops=dict(arrowstyle='<->', lw=2, color='red'))
    ax.text(stagger_time_h/2, 1.3, f'Stagger time τ = {stagger_time_h:.1f} h', 
           ha='center', va='bottom', color='red', fontweight='bold')
    
    plt.tight_layout()
    return fig

def plot_hohmann_transfer_geometry():
    """
    Plot the Hohmann transfer geometry from 1083 km to 1000 km altitude
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Constants
    moon_radius = 1737.4  # km
    r1 = moon_radius + 1083  # Initial orbit (hyperbolic perigee)
    r2 = moon_radius + 1000  # Final orbit
    
    # Create circular orbits
    theta = np.linspace(0, 2*np.pi, 100)
    
    # Moon surface
    moon_x = moon_radius * np.cos(theta)
    moon_y = moon_radius * np.sin(theta)
    
    # Orbits
    orbit1_x = r1 * np.cos(theta)
    orbit1_y = r1 * np.sin(theta)
    orbit2_x = r2 * np.cos(theta)
    orbit2_y = r2 * np.sin(theta)
    
    # Transfer ellipse
    a_transfer = (r1 + r2) / 2
    c_transfer = (r1 - r2) / 2
    b_transfer = np.sqrt(a_transfer**2 - c_transfer**2)
    
    # Ellipse parameters (transfer orbit)
    ellipse_theta = np.linspace(0, np.pi, 50)  # Half orbit
    ellipse_x = a_transfer * np.cos(ellipse_theta) + c_transfer
    ellipse_y = b_transfer * np.sin(ellipse_theta)
    
    # Plot 1: Orbital geometry
    ax1.fill(moon_x, moon_y, color='gray', alpha=0.7, label='Moon')
    ax1.plot(orbit1_x, orbit1_y, 'b-', linewidth=2, label='1083 km orbit')
    ax1.plot(orbit2_x, orbit2_y, 'g-', linewidth=2, label='1000 km orbit')
    ax1.plot(ellipse_x, ellipse_y, 'r--', linewidth=2, label='Transfer orbit')
    
    # Mark burn points
    ax1.scatter(r1, 0, color='blue', s=100, marker='o', zorder=5)
    ax1.scatter(r2, 0, color='green', s=100, marker='s', zorder=5)
    ax1.text(r1+100, 100, 'Burn 1\n(Retrograde)', ha='left', va='bottom', 
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    ax1.text(r2-100, -200, 'Burn 2\n(Prograde)', ha='right', va='top',
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    ax1.set_xlim(-3200, 3200)
    ax1.set_ylim(-3200, 3200)
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.set_title('Hohmann Transfer Geometry')
    ax1.set_xlabel('Distance (km)')
    ax1.set_ylabel('Distance (km)')
    
    # Plot 2: Velocity profile
    # Calculate velocities at different points
    mu_moon = 4902.800  # km³/s²
    
    # Circular velocities
    v1_circular = np.sqrt(mu_moon / r1)
    v2_circular = np.sqrt(mu_moon / r2)
    
    # Transfer velocities
    a_transfer = (r1 + r2) / 2
    v1_transfer = np.sqrt(mu_moon * (2/r1 - 1/a_transfer))
    v2_transfer = np.sqrt(mu_moon * (2/r2 - 1/a_transfer))
    
    # Delta-Vs
    dv1 = abs(v1_circular - v1_transfer)
    dv2 = abs(v2_circular - v2_transfer)

    positions = ['1083 km\n(Start)', '1083 km\n(Transfer)', '1000 km\n(Transfer)', '1000 km\n(Final)']
    velocities = [v1_circular, v1_transfer, v2_transfer, v2_circular]
    
    ax2.bar(positions, [v*1000 for v in velocities], color=['blue', 'red', 'red', 'green'], alpha=0.7)
    
    # Add delta-V annotations
    ax2.annotate('', xy=(0, v1_circular*1000), xytext=(1, v1_transfer*1000),
               arrowprops=dict(arrowstyle='<->', lw=2, color='purple'))
    ax2.text(0.5, (v1_circular + v1_transfer)*500, f'ΔV₁ = {dv1*1000:.1f} m/s',
           ha='center', va='bottom', color='purple', fontweight='bold')
    
    ax2.annotate('', xy=(2, v2_transfer*1000), xytext=(3, v2_circular*1000),
               arrowprops=dict(arrowstyle='<->', lw=2, color='orange'))
    ax2.text(2.5, (v2_transfer + v2_circular)*500, f'ΔV₂ = {dv2*1000:.1f} m/s',
           ha='center', va='bottom', color='orange', fontweight='bold')
    
    ax2.set_ylabel('Velocity (m/s)')
    ax2.set_title('Velocity Profile During Hohmann Transfer')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig

def plot_final_constellation():
    """
    Plot the final constellation geometry showing 60° spacing
    """
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Constants
    moon_radius = 1737.4  # km
    orbit_radius = moon_radius + 1000  # Final orbit
    
    # Moon
    theta = np.linspace(0, 2*np.pi, 100)
    moon_x = moon_radius * np.cos(theta)
    moon_y = moon_radius * np.sin(theta)
    
    # Orbit
    orbit_x = orbit_radius * np.cos(theta)
    orbit_y = orbit_radius * np.sin(theta)
    
    ax.fill(moon_x, moon_y, color='gray', alpha=0.7, label='Moon')
    ax.plot(orbit_x, orbit_y, 'k-', linewidth=2, label='1000 km orbit')
    
    # Satellite positions (60° apart)
    satellite_angles = np.arange(0, 360, 60)  # degrees
    colors = plt.cm.tab10(np.linspace(0, 1, 6))
    
    for i, angle in enumerate(satellite_angles):
        angle_rad = np.radians(angle)
        sat_x = orbit_radius * np.cos(angle_rad)
        sat_y = orbit_radius * np.sin(angle_rad)
        
        ax.scatter(sat_x, sat_y, color=colors[i], s=100, marker='o', 
                  label=f'Sat {i+1}', zorder=5)
        
        # Draw radius line
        ax.plot([0, sat_x], [0, sat_y], color=colors[i], alpha=0.5, linestyle='--')
        
        # Angle annotation
        ax.text(sat_x*1.1, sat_y*1.1, f'Sat {i+1}\n{angle}°', 
               ha='center', va='center', fontsize=8,
               bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    # Draw 60° arc between first two satellites
    arc_angles = np.linspace(0, np.radians(60), 20)
    arc_radius = orbit_radius * 0.7
    arc_x = arc_radius * np.cos(arc_angles)
    arc_y = arc_radius * np.sin(arc_angles)
    ax.plot(arc_x, arc_y, 'red', linewidth=3, label='60° spacing')
    ax.text(arc_radius*0.7, arc_radius*0.3, '60°', ha='center', va='center', 
           color='red', fontweight='bold', fontsize=12)
    
    ax.set_xlim(-3500, 3500)
    ax.set_ylim(-3500, 3500)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
    ax.set_title('Final Constellation Geometry\n6 Satellites with 60° Spacing')
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Distance (km)')
    
    plt.tight_layout()
    return fig

def plot_delta_v_budget(results):
    """
    Create a pie chart showing delta-V budget breakdown
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Extract data
    insertion_dv = results['insertion_data']['total_insertion_delta_v'] * 1000  # m/s
    hohmann_dv_total = results['deployment_schedule'][0]['total_delta_v_ms'] * 6  # All satellites
    
    # Pie chart 1: Total mission delta-V breakdown
    labels1 = ['Hyperbolic Insertion\n& Inclination Change', 'Hohmann Transfers\n(All 6 Satellites)']
    sizes1 = [insertion_dv, hohmann_dv_total]
    colors1 = ['lightcoral', 'skyblue']
    
    ax1.pie(sizes1, labels=labels1, colors=colors1, autopct='%1.1f%%', startangle=90)
    ax1.set_title(f'Total Mission Delta-V Budget\nTotal: {insertion_dv + hohmann_dv_total:.0f} m/s')
    
    # Pie chart 2: Per-satellite breakdown
    insertion_per_sat = results['insertion_data']['delta_v_insertion'] * 1000
    inclination_per_sat = results['insertion_data']['delta_v_inclination'] * 1000
    hohmann_per_sat = results['deployment_schedule'][0]['total_delta_v_ms']
    
    labels2 = ['Orbit Insertion', 'Inclination Change', 'Hohmann Transfer']
    sizes2 = [insertion_per_sat, inclination_per_sat, hohmann_per_sat]
    colors2 = ['lightgreen', 'orange', 'lightyellow']
    
    ax2.pie(sizes2, labels=labels2, colors=colors2, autopct='%1.1f%%', startangle=90)
    ax2.set_title(f'Per-Satellite Delta-V Budget\nTotal: {sum(sizes2):.0f} m/s')
    
    plt.tight_layout()
    return fig

def create_all_visualizations():
    """
    Create all visualizations for the deployment system
    """
    print("Creating visualizations for Moon Constellation Deployment...")
    
    # Run the deployment analysis
    deployer = MoonConstellationDeployer()
    results = deployer.run_deployment_analysis()
    
    print("\nGenerating plots...")
    
    # Create plots
    fig1 = plot_deployment_timeline(results['deployment_schedule'], 
                                  deployer.t_transfer_hours, 
                                  deployer.tau_hours)
    
    fig2 = plot_hohmann_transfer_geometry()
    
    fig3 = plot_final_constellation()
    
    fig4 = plot_delta_v_budget(results)
    
    # Save plots
    figs = [fig1, fig2, fig3, fig4]
    names = ['deployment_timeline', 'hohmann_geometry', 'final_constellation', 'delta_v_budget']
    
    for fig, name in zip(figs, names):
        filename = f'moon_constellation_{name}.png'
        fig.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Saved: {filename}")
    
    # Show all plots
    plt.show()
    
    return results

if __name__ == "__main__":
    results = create_all_visualizations()