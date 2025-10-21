"""
Satellite Deployment System for Moon Communications Constellation
Starting from Polar Circular Orbit at 1084 km

This system implements the deployment strategy:
1. Start from existing polar circular orbit at 1084 km altitude
2. Deploy 6 satellites with 60° spacing in the 1084km orbit
3. Transfer all satellites down to 1000km final orbit using Hohmann transfers

Based on lunar constants:
- Moon radius R = 1737.4 km
- Moon gravitational parameter μ = 4902.800 km³/s²
"""

import math
import sys
import os

# Add Transfer_Simulation to path
transfer_sim_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'Transfer_Simulation')
sys.path.insert(0, transfer_sim_path)

# Lunar constants (from the images)
MOON_RADIUS_KM = 1737.4  # km
MU_MOON_KMS = 4902.800   # km³/s²

# Mission parameters
R1 = MOON_RADIUS_KM + 1084  # Initial polar circular orbit (1084 km altitude)
R2 = MOON_RADIUS_KM + 1000  # Final circular orbit (1000 km altitude)
TARGET_SPACING_DEG = 60     # degrees between satellites
NUM_SATELLITES = 6          # Total number of satellites
TOTAL_SATELLITES = NUM_SATELLITES  # 6 total satellites

class MoonConstellationDeployer:
    """
    Handles deployment of 6 satellites starting from polar circular orbit at 1084km:
    1. Start from existing polar circular orbit at 1084km altitude (90° inclination)
    2. Deploy 6 satellites with 60° spacing in the 1084km orbit
    3. Transfer all satellites down to 1000km final orbit using Hohmann transfers
    """
    
    def __init__(self):
        self.moon_radius = MOON_RADIUS_KM
        self.mu_moon = MU_MOON_KMS
        self.r1 = R1  # Initial orbit radius (1084 km altitude - polar circular)
        self.r2 = R2  # Final orbit radius (1000 km altitude)
        self.target_spacing_deg = TARGET_SPACING_DEG
        self.num_satellites = NUM_SATELLITES
        self.total_satellites = TOTAL_SATELLITES
        
        # Calculate orbital parameters
        self.calculate_orbital_parameters()
        self.calculate_transfer_parameters()
        self.calculate_staggering_time()
    
    def calculate_orbital_parameters(self):
        """Calculate circular orbital parameters for both orbits"""
        # Initial orbit (1083 km altitude)
        self.n1 = math.sqrt(self.mu_moon / self.r1**3)  # Mean motion rad/s
        self.v1_circular = math.sqrt(self.mu_moon / self.r1)  # Circular velocity km/s
        self.T1 = 2 * math.pi / self.n1  # Orbital period seconds
        self.T1_hours = self.T1 / 3600  # Period in hours
        
        # Final orbit (1000 km altitude)
        self.n2 = math.sqrt(self.mu_moon / self.r2**3)  # Mean motion rad/s
        self.v2_circular = math.sqrt(self.mu_moon / self.r2)  # Circular velocity km/s
        self.T2 = 2 * math.pi / self.n2  # Orbital period seconds
        self.T2_hours = self.T2 / 3600  # Period in hours
        
        print(f"Orbital Parameters:")
        print(f"  Initial orbit (1084 km): n1 = {self.n1:.4e} rad/s, T1 = {self.T1_hours:.2f} h")
        print(f"  Final orbit (1000 km):   n2 = {self.n2:.4e} rad/s, T2 = {self.T2_hours:.2f} h")
    
    def calculate_transfer_parameters(self):
        """Calculate Hohmann transfer parameters"""
        # Transfer ellipse parameters
        self.a_transfer = (self.r1 + self.r2) / 2  # Semi-major axis
        
        # Velocities for Hohmann transfer
        # At 1084 km (apogee of transfer orbit)
        self.v_transfer_r1 = math.sqrt(self.mu_moon * (2/self.r1 - 1/self.a_transfer))
        
        # At 1000 km (perigee of transfer orbit)  
        self.v_transfer_r2 = math.sqrt(self.mu_moon * (2/self.r2 - 1/self.a_transfer))
        
        # Delta-V calculations
        self.delta_v1 = abs(self.v1_circular - self.v_transfer_r1)  # First burn (retrograde at 1084 km)
        self.delta_v2 = abs(self.v2_circular - self.v_transfer_r2)  # Second burn (prograde at 1000 km)
        self.total_delta_v = self.delta_v1 + self.delta_v2
        
        # Transfer time (half period of transfer ellipse)
        self.t_transfer = math.pi * math.sqrt(self.a_transfer**3 / self.mu_moon)
        self.t_transfer_hours = self.t_transfer / 3600
        
        print(f"\nHohmann Transfer Parameters:")
        print(f"  Transfer semi-major axis: {self.a_transfer:.1f} km")
        print(f"  First burn (retrograde at 1084 km): Delta-v1 = {self.delta_v1*1000:.2f} m/s")
        print(f"  Second burn (prograde at 1000 km): Delta-v2 = {self.delta_v2*1000:.2f} m/s")
        print(f"  Total per satellite: Delta-v = {self.total_delta_v*1000:.2f} m/s")
        print(f"  Transfer time: {self.t_transfer_hours:.3f} h ({self.t_transfer:.1f} s)")
    
    def calculate_staggering_time(self):
        """Calculate the staggering time τ for 60° spacing in final orbit"""
        # The staggering time is based on differential motion between initial and final orbits
        # Δθₖ = (n₂ - n₁)(k-1)τ
        # For 60° spacing: Δθ = π/3 radians
        # τ = π/3 / (n₂ - n₁)
        
        spacing_rad = math.radians(self.target_spacing_deg)
        self.tau = spacing_rad / (self.n2 - self.n1)  # Correct formula restored
        self.tau_hours = self.tau / 3600
        
        print(f"\nStaggering Parameters:")
        print(f"  Target spacing: {self.target_spacing_deg} deg = {spacing_rad:.4f} rad")
        print(f"  Mean motion difference: n2 - n1 = {(self.n2 - self.n1):.4e} rad/s")
        print(f"  Staggering time tau = {self.tau:.1f} s = {self.tau_hours:.3f} h")
    
    def get_initial_orbit_data(self):
        """
        Get data for the initial polar circular orbit at 1084 km
        
        Returns:
        - Dictionary with initial orbit data
        """
        # Calculate velocity for circular orbit at 1084 km altitude
        v_circular = math.sqrt(self.mu_moon / self.r1)
        
        # Orbital period
        T_circular = 2 * math.pi * math.sqrt(self.r1**3 / self.mu_moon)
        
        initial_orbit_data = {
            'altitude_km': 1084,
            'radius_km': self.r1,
            'velocity_kms': v_circular,
            'period_s': T_circular,
            'period_h': T_circular / 3600,
            'inclination_deg': 90.0,
            'orbit_type': 'Polar Circular',
            'mean_motion_rad_s': self.n1
        }
        
        return initial_orbit_data
    
    def calculate_raan_change_maneuver(self, raan_change_deg):
        """
        Calculate delta-V required to change RAAN without changing inclination.
        This maneuver must be performed at the polar crossing points (90° from ascending node).
        
        Parameters:
        - raan_change_deg: RAAN change in degrees
        
        Returns:
        - Dictionary with RAAN change maneuver data
        """
        # For a polar orbit (i = 90°), RAAN change is performed at the polar crossings
        # At these points, the spacecraft is moving purely in the orbital plane
        
        # Convert RAAN change to radians
        raan_change_rad = math.radians(raan_change_deg)
        
        # Velocity at 1084km circular orbit
        v_circular = self.v1_circular  # km/s
        
        # For RAAN change at polar crossing in a polar orbit:
        # The velocity vector is horizontal (tangent to surface)
        # Delta-V = 2 * v * sin(ΔRAAN/2) 
        delta_v_raan = 2 * v_circular * math.sin(raan_change_rad / 2)
        
        # Time to reach polar crossing from ascending node
        # For a polar orbit, polar crossings are at 90° and 270° true anomaly
        time_to_polar_crossing = (math.pi / 2) / self.n1  # Quarter orbit
        time_to_polar_crossing_hours = time_to_polar_crossing / 3600
        
        raan_data = {
            'raan_change_deg': raan_change_deg,
            'raan_change_rad': raan_change_rad,
            'delta_v_raan_ms': delta_v_raan * 1000,
            'delta_v_raan_kms': delta_v_raan,
            'maneuver_location': 'Polar crossing (90° from ascending node)',
            'orbital_velocity_kms': v_circular,
            'time_to_maneuver_s': time_to_polar_crossing,
            'time_to_maneuver_h': time_to_polar_crossing_hours,
            'altitude_km': 1084,
            'inclination_deg': 90
        }
        
        print(f"\nRAAN CHANGE MANEUVER ANALYSIS:")
        print(f"  RAAN change required: {raan_change_deg:.2f} deg = {raan_change_rad:.4f} rad")
        print(f"  Maneuver location: Polar crossing (±90° latitude)")
        print(f"  Orbital velocity: {v_circular:.3f} km/s")
        print(f"  Delta-V required: {delta_v_raan*1000:.2f} m/s")
        print(f"  Time to polar crossing: {time_to_polar_crossing_hours:.3f} h ({time_to_polar_crossing:.1f} s)")
        print(f"  Maneuver efficiency: Optimal (performed at polar crossing)")
        
        return raan_data
    
    def calculate_multi_plane_deployment(self, num_planes=7):
        """
        Calculate deployment strategy for multiple orbital planes using RAAN changes.
        
        Parameters:
        - num_planes: Number of orbital planes (default 7)
        
        Returns:
        - Dictionary with multi-plane deployment analysis
        """
        raan_separation_deg = 360.0 / num_planes
        
        print(f"\nMULTI-PLANE DEPLOYMENT ANALYSIS:")
        print(f"  Number of orbital planes: {num_planes}")
        print(f"  RAAN separation: {raan_separation_deg:.2f} deg")
        
        planes_data = []
        total_raan_delta_v = 0
        
        for plane_id in range(num_planes):
            if plane_id == 0:
                # First plane uses the original RAAN (no maneuver needed)
                raan_change_deg = 0
                delta_v_raan = 0
            else:
                # Subsequent planes need RAAN change
                raan_change_deg = raan_separation_deg
                raan_maneuver = self.calculate_raan_change_maneuver(raan_change_deg)
                delta_v_raan = raan_maneuver['delta_v_raan_ms']
                total_raan_delta_v += delta_v_raan
            
            plane_raan_deg = plane_id * raan_separation_deg
            
            plane_data = {
                'plane_id': plane_id + 1,
                'plane_raan_deg': plane_raan_deg,
                'raan_change_from_previous_deg': raan_change_deg,
                'delta_v_raan_ms': delta_v_raan,
                'satellites_per_plane': 1  # Could be modified for multiple sats per plane
            }
            
            planes_data.append(plane_data)
            
            print(f"  Plane {plane_id + 1}: RAAN = {plane_raan_deg:.1f}°, "
                  f"ΔV = {delta_v_raan:.1f} m/s")
        
        print(f"\nMULTI-PLANE SUMMARY:")
        print(f"  Total RAAN change Delta-V: {total_raan_delta_v:.1f} m/s")
        print(f"  Average Delta-V per plane: {total_raan_delta_v/(num_planes-1):.1f} m/s")
        print(f"  Coverage: {num_planes} equally-spaced orbital planes")
        
        return {
            'num_planes': num_planes,
            'raan_separation_deg': raan_separation_deg,
            'planes_data': planes_data,
            'total_raan_delta_v_ms': total_raan_delta_v,
            'average_raan_delta_v_ms': total_raan_delta_v/(num_planes-1) if num_planes > 1 else 0
        }
    
    def create_deployment_schedule(self):
        """
        Create deployment schedule for 6 satellites:
        1. Deploy all 6 satellites together at the same time in 1084km orbit
        2. Use staggered Hohmann transfers to achieve 60° spacing in final 1000km orbit
        """
        deployment_schedule = {
            'phase1_deployment': {
                'description': 'Deploy all satellites together in 1084km orbit',
                'time_s': 0,
                'time_h': 0,
                'satellites': []
            },
            'phase2_transfer': {
                'description': 'Staggered transfers to achieve 60° spacing in 1000km orbit',
                'satellites': []
            }
        }
        
        # Phase 1: Deploy all satellites together at t=0
        for sat_id in range(1, self.num_satellites + 1):
            satellite_deployment = {
                'satellite_id': sat_id,
                'deployment_time_s': 0,  # All deployed at t=0
                'deployment_time_h': 0,  # All deployed at t=0
                'initial_angular_position_deg': 0,  # All start at same position
                'orbit_altitude_km': 1084,
                'orbit_inclination_deg': 90
            }
            
            deployment_schedule['phase1_deployment']['satellites'].append(satellite_deployment)
        
        # Phase 2: Staggered transfers to achieve spacing
        # Satellites are transferred with time intervals τ to achieve 60° spacing in final orbit
        stabilization_time = self.T1  # One full orbit for stabilization after deployment
        
        for sat_id in range(1, self.num_satellites + 1):
            # Each satellite starts transfer at staggered intervals
            # The staggering creates the desired angular separation in the final orbit
            transfer_start_time = stabilization_time + (sat_id - 1) * self.tau
            
            # Transfer burn times
            first_burn_time = transfer_start_time
            second_burn_time = transfer_start_time + self.t_transfer
            
            # Final angular position due to staggered transfers
            final_angular_position_deg = (sat_id - 1) * self.target_spacing_deg
            
            satellite_transfer = {
                'satellite_id': sat_id,
                'transfer_start_time_s': transfer_start_time,
                'transfer_start_time_h': transfer_start_time / 3600,
                'first_burn_time_s': first_burn_time,
                'first_burn_time_h': first_burn_time / 3600,
                'second_burn_time_s': second_burn_time,
                'second_burn_time_h': second_burn_time / 3600,
                'delta_v1_ms': self.delta_v1 * 1000,
                'delta_v2_ms': self.delta_v2 * 1000,
                'total_transfer_delta_v_ms': self.total_delta_v * 1000,
                'final_altitude_km': 1000,
                'final_inclination_deg': 90,
                'final_angular_position_deg': final_angular_position_deg
            }
            
            deployment_schedule['phase2_transfer']['satellites'].append(satellite_transfer)
        
        return deployment_schedule
    
    def display_mission_analysis(self, deployment_schedule):
        """Display comprehensive mission analysis for constellation deployment"""
        print("\n" + "=" * 80)
        print("MOON COMMUNICATIONS CONSTELLATION DEPLOYMENT ANALYSIS")
        print("=" * 80)
        
        print(f"\nMISSION OVERVIEW:")
        print(f"   Constellation: {self.total_satellites} satellites in single orbital plane")
        print(f"   Deployment orbit: 1084 km circular, 90 deg inclination")
        print(f"   Final orbits: 1000 km circular, 90 deg inclination")
        print(f"   Spacing strategy: {self.target_spacing_deg} deg between satellites")
        print(f"   Deployment method: Staggered deployment with Hohmann transfers")
        
        print(f"\nPHASE 1: SATELLITE DEPLOYMENT IN 1084KM ORBIT")
        print(f"   Deployment orbit: 1084 km altitude, 90 deg inclination")
        print(f"   Number of satellites: {self.num_satellites}")
        print(f"   Deployment method: All satellites deployed together at t=0")
        print(f"   Initial spacing: 0° (all satellites co-located initially)")
        print(f"   Orbit period: {self.T1_hours:.2f} h")
        
        # Show deployment timeline
        print(f"\nDEPLOYMENT TIMELINE (1084km orbit):")
        print(f"   All {self.num_satellites} satellites deployed simultaneously at t=0")
        print(f"   Initial angular position: 0° for all satellites")
        print(f"   Stabilization period: {self.T1_hours:.2f} hours (1 orbit)")
        
        print(f"\nPHASE 2: STAGGERED HOHMANN TRANSFERS TO 1000KM ORBIT")
        print(f"   Transfer method: Hohmann transfers (1084 km -> 1000 km)")
        print(f"   Staggering strategy: Transfers spaced by τ = {self.tau_hours:.3f} h ({self.tau:.1f} s)")
        print(f"   Transfer time: {self.t_transfer_hours:.3f} h per satellite")
        print(f"   First burn Delta-v: {self.delta_v1*1000:.2f} m/s (retrograde at 1084 km)")
        print(f"   Second burn Delta-v: {self.delta_v2*1000:.2f} m/s (prograde at 1000 km)")
        print(f"   Total transfer Delta-v: {self.total_delta_v*1000:.2f} m/s per satellite")
        print(f"   Spacing mechanism: Staggered transfers create 60° spacing in final orbit")
        
        # Show transfer timeline
        print(f"\nSTAGGERED TRANSFER TIMELINE (1084km -> 1000km):")
        print(f"   {'Sat ID':<8} {'Start Time':<12} {'1st Burn':<12} {'2nd Burn':<12} {'Final Pos':<12}")
        print(f"   {'#':<8} {'(hours)':<12} {'(hours)':<12} {'(hours)':<12} {'(degrees)':<12}")
        print(f"   {'-'*8} {'-'*12} {'-'*12} {'-'*12} {'-'*12}")
        
        total_mission_time = 0
        for sat in deployment_schedule['phase2_transfer']['satellites']:
            total_mission_time = max(total_mission_time, sat['second_burn_time_h'])
            print(f"   {sat['satellite_id']:<8} {sat['transfer_start_time_h']:<12.3f} "
                  f"{sat['first_burn_time_h']:<12.3f} {sat['second_burn_time_h']:<12.3f} "
                  f"{sat['final_angular_position_deg']:<12.1f}")
        
        # Delta-V budget (only transfer maneuvers needed)
        total_transfer_delta_v = self.num_satellites * self.total_delta_v * 1000
        total_mission_delta_v = total_transfer_delta_v  # No insertion maneuvers needed
        
        print(f"\nDELTA-V BUDGET:")
        print(f"   Initial orbit setup: 0.0 m/s (already in polar circular orbit)")
        print(f"   Hohmann transfers ({self.num_satellites} satellites): {total_transfer_delta_v:.1f} m/s")
        print(f"   Total mission Delta-v: {total_mission_delta_v:.1f} m/s")
        
        print(f"\nMISSION TIMELINE:")
        phase1_completion = self.T1_hours  # Just stabilization time
        print(f"   Phase 1 completion: {phase1_completion:.2f} hours (deployment + stabilization)")
        print(f"   Phase 2 completion: {total_mission_time:.2f} hours (all satellites in final orbit)")
        print(f"   Total mission duration: {total_mission_time:.2f} hours")
        
        # Calculate phase durations for analysis
        deployment_duration = phase1_completion  # Just stabilization
        staggered_transfer_duration = total_mission_time - phase1_completion
        
        print(f"\nPHASE DURATION ANALYSIS:")
        print(f"   Deployment + stabilization: {deployment_duration:.2f} hours ({deployment_duration/24:.1f} days)")
        print(f"   Staggered transfer phase: {staggered_transfer_duration:.2f} hours ({staggered_transfer_duration/24:.1f} days)")
        print(f"   Transfer staggering span: {(self.num_satellites-1) * self.tau_hours:.2f} hours")
        print(f"   Individual transfer time: {self.t_transfer_hours:.2f} hours per satellite")
        
        print(f"\nFINAL CONSTELLATION:")
        print(f"   Total satellites: {self.total_satellites}")
        print(f"   Orbit altitude: 1000 km")
        print(f"   Orbit inclination: 90 deg (polar)")
        print(f"   Orbital planes: 1")
        print(f"   Satellites per plane: {self.num_satellites} (separated by {self.target_spacing_deg} deg)")
        print(f"   Coverage: Polar lunar coverage")
        print(f"   Final orbit period: {self.T2_hours:.2f} hours")
        print(f"   Constellation revisit time: ~{self.T2_hours/self.num_satellites:.2f} hours")
        
        return {
            'deployment_schedule': deployment_schedule,
            'mission_summary': {
                'total_mission_time_h': total_mission_time,
                'deployment_phase_h': deployment_duration,
                'staggered_transfer_phase_h': staggered_transfer_duration,
                'total_mission_delta_v_ms': total_mission_delta_v,
                'total_satellites': self.total_satellites,
                'deployment_orbit_km': 1084,
                'final_orbit_km': 1000
            }
        }
    
    def run_deployment_analysis(self):
        """
        Run complete constellation deployment analysis starting from polar circular orbit at 1084km
        """
        print("MOON CONSTELLATION DEPLOYER")
        print("Deployment from 1084km Polar Circular Orbit")
        print("=" * 50)
        
        # Create deployment schedule
        deployment_schedule = self.create_deployment_schedule()
        
        # Display complete analysis
        results = self.display_mission_analysis(deployment_schedule)
        
        return results


def run_standard_deployment():
    """
    Run standard deployment analysis starting from polar circular orbit at 1084km
    """
    print("Starting constellation deployment from existing polar circular orbit...")
    
    deployer = MoonConstellationDeployer()
    results = deployer.run_deployment_analysis()
    
    return results


def test_raan_change_analysis():
    """Test function for RAAN change analysis"""
    print("RAAN CHANGE ANALYSIS TEST")
    print("=" * 50)
    
    deployer = MoonConstellationDeployer()
    
    # Calculate RAAN change for 360/7 degrees
    raan_change_deg = 360.0 / 7
    raan_data = deployer.calculate_raan_change_maneuver(raan_change_deg)
    
    # Calculate multi-plane deployment
    multi_plane_data = deployer.calculate_multi_plane_deployment(num_planes=7)
    
    return raan_data, multi_plane_data


def main():
    """Main function"""
    print("Select analysis mode:")
    print("1. Standard constellation deployment (from 1084km polar orbit)")
    print("2. RAAN change analysis (multi-plane deployment)")
    
    try:
        choice = input("Enter choice (1 or 2): ").strip()
        
        if choice == '2':
            results = test_raan_change_analysis()
        else:
            results = run_standard_deployment()
            
        return results
        
    except KeyboardInterrupt:
        print("\nOperation cancelled by user")
        return None
    except Exception as e:
        print(f"Error: {e}")
        # Fallback to standard deployment
        results = run_standard_deployment()
        return results


if __name__ == "__main__":
    results = main()
