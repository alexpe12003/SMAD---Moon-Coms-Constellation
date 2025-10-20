"""
Satellite Deployment System for Moon Communications Constellation
Using Hohmann Transfer Strategy with Staggered Burns

This system implements the deployment strategy shown in the reference images:
1. Circularize hyperbolic orbit at perigee and change inclination to 90°
2. Use Hohmann transfers from 1100km to 1000km circular orbits
3. Stagger burns with time gap τ to achieve 60° spacing between satellites

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
R1 = MOON_RADIUS_KM + 1764  # Initial circular orbit at hyperbolic perigee (1764 km altitude)
R2 = MOON_RADIUS_KM + 1000  # Final circular orbit (1000 km altitude)
TARGET_SPACING_DEG = 60     # degrees between satellites within each carrier
NUM_CARRIERS = 7            # Number of carriers
NUM_SATELLITES_PER_CARRIER = 6  # Satellites per carrier
TOTAL_SATELLITES = NUM_CARRIERS * NUM_SATELLITES_PER_CARRIER  # 42 total satellites

class MoonConstellationDeployer:
    """
    Handles deployment of 42 satellites (7 carriers × 6 satellites each) using:
    1. Staggered inclination changes for carriers (different RAAN)
    2. Staggered Hohmann transfers within each carrier
    """
    
    def __init__(self):
        self.moon_radius = MOON_RADIUS_KM
        self.mu_moon = MU_MOON_KMS
        self.r1 = R1  # Initial orbit radius (1764 km altitude - hyperbolic perigee)
        self.r2 = R2  # Final orbit radius (1000 km altitude)
        self.target_spacing_deg = TARGET_SPACING_DEG
        self.num_carriers = NUM_CARRIERS
        self.num_satellites_per_carrier = NUM_SATELLITES_PER_CARRIER
        self.total_satellites = TOTAL_SATELLITES
        
        # Calculate orbital parameters
        self.calculate_orbital_parameters()
        self.calculate_transfer_parameters()
        self.calculate_staggering_time()
        self.calculate_carrier_staggering_time()
    
    def calculate_orbital_parameters(self):
        """Calculate circular orbital parameters for both orbits"""
        # Initial orbit (1100 km altitude)
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
        print(f"  Initial orbit (1764 km): n1 = {self.n1:.4e} rad/s, T1 = {self.T1_hours:.2f} h")
        print(f"  Final orbit (1000 km):   n2 = {self.n2:.4e} rad/s, T2 = {self.T2_hours:.2f} h")
    
    def calculate_transfer_parameters(self):
        """Calculate Hohmann transfer parameters"""
        # Transfer ellipse parameters
        self.a_transfer = (self.r1 + self.r2) / 2  # Semi-major axis
        
        # Velocities for Hohmann transfer
        # At 1764 km (apogee of transfer orbit)
        self.v_transfer_r1 = math.sqrt(self.mu_moon * (2/self.r1 - 1/self.a_transfer))
        
        # At 1000 km (perigee of transfer orbit)  
        self.v_transfer_r2 = math.sqrt(self.mu_moon * (2/self.r2 - 1/self.a_transfer))
        
        # Delta-V calculations
        self.delta_v1 = abs(self.v1_circular - self.v_transfer_r1)  # First burn (retrograde at 1764 km)
        self.delta_v2 = abs(self.v2_circular - self.v_transfer_r2)  # Second burn (prograde at 1000 km)
        self.total_delta_v = self.delta_v1 + self.delta_v2
        
        # Transfer time (half period of transfer ellipse)
        self.t_transfer = math.pi * math.sqrt(self.a_transfer**3 / self.mu_moon)
        self.t_transfer_hours = self.t_transfer / 3600
        
        print(f"\nHohmann Transfer Parameters:")
        print(f"  Transfer semi-major axis: {self.a_transfer:.1f} km")
        print(f"  First burn (retrograde at 1764 km): Delta-v1 = {self.delta_v1*1000:.2f} m/s")
        print(f"  Second burn (prograde at 1000 km): Delta-v2 = {self.delta_v2*1000:.2f} m/s")
        print(f"  Total per satellite: Delta-v = {self.total_delta_v*1000:.2f} m/s")
        print(f"  Transfer time: {self.t_transfer_hours:.3f} h ({self.t_transfer:.1f} s)")
    
    def calculate_staggering_time(self):
        """Calculate the staggering time τ for 60° spacing"""
        # From the images: Δθₖ = (n₁ - n₂)(k-1)τ
        # For 60° spacing: Δθ = π/3 radians
        # τ = π/3 / (n₂ - n₁)
        
        spacing_rad = math.radians(self.target_spacing_deg)
        self.tau = spacing_rad / (self.n2 - self.n1)
        self.tau_hours = self.tau / 3600
        
        print(f"\nStaggering Parameters:")
        print(f"  Target spacing: {self.target_spacing_deg} deg = {spacing_rad:.4f} rad")
        print(f"  Mean motion difference: n2 - n1 = {(self.n2 - self.n1):.4e} rad/s")
        print(f"  Staggering time tau = {self.tau:.1f} s = {self.tau_hours:.3f} h")
    
    def calculate_carrier_staggering_time(self):
        """Calculate the staggering time for carrier inclination changes"""
        # Each carrier waits T1/7 before inclination change to get different RAAN
        self.carrier_stagger_time = self.T1 / self.num_carriers  # seconds
        self.carrier_stagger_hours = self.carrier_stagger_time / 3600  # hours
        
        print(f"\nCarrier Staggering Parameters:")
        print(f"  Number of carriers: {self.num_carriers}")
        print(f"  Initial orbit period T1: {self.T1_hours:.3f} h")
        print(f"  Carrier stagger time: T1/{self.num_carriers} = {self.carrier_stagger_hours:.3f} h ({self.carrier_stagger_time:.1f} s)")
        print(f"  RAAN separation: {360/self.num_carriers:.1f} deg per carrier")
    
    def calculate_hyperbolic_insertion_maneuver(self, hyperbolic_data):
        """
        Calculate maneuver to insert from hyperbolic orbit into circular orbit at perigee (1764 km)
        and change inclination to 90°
        
        Parameters:
        - hyperbolic_data: Dictionary with hyperbolic orbit parameters
        
        Returns:
        - Dictionary with insertion maneuver data
        """
        # From Transfer_Simulation results, we get hyperbolic velocity at perigee
        v_hyp_perigee = hyperbolic_data.get('v_hyp_perigee', 2.077)  # km/s from previous run
        perigee_altitude = hyperbolic_data.get('perigee_altitude', 1764)  # km
        
        # Calculate velocity needed for circular orbit at perigee altitude (1764 km)
        v_circular_perigee = math.sqrt(self.mu_moon / self.r1)
        
        # Delta-V for orbit insertion (velocity change magnitude)
        delta_v_insertion = abs(v_hyp_perigee - v_circular_perigee)
        
        # Delta-V for inclination change from 28.6° to 90°
        # Δv = 2*v*sin(Δi/2) where Δi = |90° - 28.6°| = 61.4°
        initial_inclination_deg = 28.6
        target_inclination_deg = 90.0
        inclination_change_deg = abs(target_inclination_deg - initial_inclination_deg)
        delta_v_inclination = 2 * v_circular_perigee * math.sin(math.radians(inclination_change_deg)/2)
        
        # Total insertion delta-v
        total_insertion_delta_v = delta_v_insertion + delta_v_inclination
        
        insertion_data = {
            'v_hyperbolic': v_hyp_perigee,
            'v_circular_perigee': v_circular_perigee,
            'delta_v_insertion': delta_v_insertion,
            'delta_v_inclination': delta_v_inclination,
            'total_insertion_delta_v': total_insertion_delta_v,
            'initial_inclination_deg': initial_inclination_deg,
            'final_inclination_deg': target_inclination_deg,
            'inclination_change_deg': inclination_change_deg,
            'maneuver_altitude': perigee_altitude  # Circularize at perigee altitude
        }
        
        return insertion_data
    
    def create_deployment_schedule(self):
        """Create the deployment schedule for all 7 carriers with 6 satellites each"""
        all_schedules = []
        
        for carrier_id in range(1, self.num_carriers + 1):
            # Time when this carrier performs inclination change
            carrier_inclination_time = (carrier_id - 1) * self.carrier_stagger_time
            carrier_inclination_hours = carrier_inclination_time / 3600
            
            # RAAN for this carrier (equally spaced)
            carrier_raan_deg = (carrier_id - 1) * (360 / self.num_carriers)
            
            carrier_schedule = {
                'carrier_id': carrier_id,
                'inclination_change_time_s': carrier_inclination_time,
                'inclination_change_time_h': carrier_inclination_hours,
                'raan_deg': carrier_raan_deg,
                'satellites': []
            }
            
            # Create satellite deployment schedule within this carrier
            for sat_id in range(1, self.num_satellites_per_carrier + 1):
                # Time for satellite deployment (relative to carrier's inclination change)
                sat_deployment_time = carrier_inclination_time + (sat_id - 1) * self.tau
                sat_deployment_hours = sat_deployment_time / 3600
                
                # Satellite burns
                first_burn_time = sat_deployment_time
                second_burn_time = sat_deployment_time + self.t_transfer
                
                # Angular separation within carrier
                sat_separation_deg = (sat_id - 1) * self.target_spacing_deg
                
                # Global satellite ID
                global_sat_id = (carrier_id - 1) * self.num_satellites_per_carrier + sat_id
                
                satellite_schedule = {
                    'global_satellite_id': global_sat_id,
                    'carrier_id': carrier_id,
                    'satellite_id_in_carrier': sat_id,
                    'deployment_time_s': sat_deployment_time,
                    'deployment_time_h': sat_deployment_hours,
                    'first_burn_time_s': first_burn_time,
                    'first_burn_time_h': first_burn_time / 3600,
                    'second_burn_time_s': second_burn_time,
                    'second_burn_time_h': second_burn_time / 3600,
                    'separation_in_carrier_deg': sat_separation_deg,
                    'raan_deg': carrier_raan_deg,
                    'delta_v1_ms': self.delta_v1 * 1000,
                    'delta_v2_ms': self.delta_v2 * 1000,
                    'total_delta_v_ms': self.total_delta_v * 1000
                }
                
                carrier_schedule['satellites'].append(satellite_schedule)
            
            all_schedules.append(carrier_schedule)
        
        return all_schedules
    
    def display_mission_analysis(self, insertion_data, deployment_schedules):
        """Display comprehensive mission analysis for multi-carrier deployment"""
        print("\n" + "=" * 80)
        print("MOON COMMUNICATIONS CONSTELLATION DEPLOYMENT ANALYSIS")
        print("=" * 80)
        
        print(f"\nMISSION OVERVIEW:")
        print(f"   Constellation: {self.total_satellites} satellites ({self.num_carriers} carriers × {self.num_satellites_per_carrier} satellites)")
        print(f"   Final orbits: 1000 km circular, 90 deg inclination")
        print(f"   Spacing strategy: {self.target_spacing_deg} deg within carriers, {360/self.num_carriers:.1f} deg RAAN separation")
        print(f"   Deployment method: Staggered carriers + Staggered Hohmann transfers")
        
        print(f"\nPHASE 1: HYPERBOLIC ORBIT INSERTION")
        print(f"   Initial state: Hyperbolic approach orbit (28.6 deg inclination)")
        print(f"   Target: 1764 km circular (at perigee)")
        print(f"   Insertion Delta-v: {insertion_data['delta_v_insertion']*1000:.1f} m/s")
        
        print(f"\nPHASE 2: CARRIER INCLINATION CHANGES")
        print(f"   Inclination change: {insertion_data['inclination_change_deg']:.1f} deg (28.6 deg -> 90 deg)")
        print(f"   Inclination Delta-v: {insertion_data['delta_v_inclination']*1000:.1f} m/s")
        print(f"   Carrier staggering: {self.carrier_stagger_hours:.3f} h between carriers")
        print(f"   RAAN separation: {360/self.num_carriers:.1f} deg per carrier")
        
        print(f"\nPHASE 3: SATELLITE DEPLOYMENT PER CARRIER")
        print(f"   Method: Staggered Hohmann transfers (1764 km -> 1000 km)")
        print(f"   Satellite staggering: {self.tau_hours:.3f} h between satellites within carrier")
        print(f"   Transfer time: {self.t_transfer_hours:.3f} h per satellite")
        
        # Show carrier summary
        print(f"\nCARRIER DEPLOYMENT SCHEDULE:")
        print(f"   {'Carrier':<8} {'RAAN':<8} {'Inc.Change':<12} {'First Sat':<12} {'Last Sat':<12}")
        print(f"   {'ID':<8} {'(deg)':<8} {'Time (h)':<12} {'Deploy (h)':<12} {'Deploy (h)':<12}")
        print(f"   {'-'*8} {'-'*8} {'-'*12} {'-'*12} {'-'*12}")
        
        total_mission_time = 0
        for carrier in deployment_schedules:
            first_sat_time = carrier['satellites'][0]['deployment_time_h']
            last_sat_time = carrier['satellites'][-1]['second_burn_time_h']
            total_mission_time = max(total_mission_time, last_sat_time)
            
            print(f"   {carrier['carrier_id']:<8} {carrier['raan_deg']:<8.1f} "
                  f"{carrier['inclination_change_time_h']:<12.3f} "
                  f"{first_sat_time:<12.3f} {last_sat_time:<12.3f}")
        
        # Delta-V budget
        hohmann_total = sum(len(carrier['satellites']) for carrier in deployment_schedules) * self.total_delta_v * 1000
        
        print(f"\nDELTA-V BUDGET:")
        print(f"   Hyperbolic insertion (all carriers): {insertion_data['delta_v_insertion']*1000:.1f} m/s")
        print(f"   Inclination changes ({self.num_carriers} carriers): {insertion_data['delta_v_inclination']*1000*self.num_carriers:.1f} m/s")
        print(f"   Hohmann transfers ({self.total_satellites} satellites): {hohmann_total:.1f} m/s")
        print(f"   Total constellation Delta-v: {(insertion_data['delta_v_insertion']*1000 + insertion_data['delta_v_inclination']*1000*self.num_carriers + hohmann_total):.1f} m/s")
        
        print(f"\nMISSION TIMELINE:")
        print(f"   Total deployment duration: {total_mission_time:.2f} hours")
        print(f"   First satellite operational: {deployment_schedules[0]['satellites'][0]['second_burn_time_h']:.3f} hours")
        print(f"   Last satellite operational: {total_mission_time:.2f} hours")
        
        print(f"\nFINAL CONSTELLATION:")
        print(f"   Total satellites: {self.total_satellites}")
        print(f"   Orbit altitude: 1000 km")
        print(f"   Orbit inclination: 90 deg (polar)")
        print(f"   Orbital planes: {self.num_carriers} (RAAN separated by {360/self.num_carriers:.1f} deg)")
        print(f"   Satellites per plane: {self.num_satellites_per_carrier} (separated by {self.target_spacing_deg} deg)")
        print(f"   Coverage: Enhanced global lunar coverage")
        print(f"   Orbital period: {self.T2_hours:.2f} hours")
        print(f"   Constellation revisit time: ~{self.T2_hours/(self.total_satellites/6):.2f} hours")
        
        return {
            'insertion_data': insertion_data,
            'deployment_schedules': deployment_schedules,
            'mission_summary': {
                'total_mission_time_h': total_mission_time,
                'total_constellation_delta_v_ms': (insertion_data['delta_v_insertion']*1000 + 
                                                  insertion_data['delta_v_inclination']*1000*self.num_carriers + 
                                                  hohmann_total),
                'total_satellites': self.total_satellites,
                'num_carriers': self.num_carriers
            }
        }
    
    def run_deployment_analysis(self, hyperbolic_orbit_data=None):
        """
        Run complete multi-carrier deployment analysis
        
        Parameters:
        - hyperbolic_orbit_data: Optional dictionary with hyperbolic orbit parameters
        """
        print("MOON CONSTELLATION DEPLOYER")
        print("Multi-Carrier Staggered Deployment Strategy")
        print("=" * 50)
        
        # Default hyperbolic data if not provided (from Transfer_Simulation output)
        if hyperbolic_orbit_data is None:
            hyperbolic_orbit_data = {
                'v_hyp_perigee': 2.077,  # km/s
                'perigee_altitude': 1764,  # km
                'eccentricity': 2.081
            }
        
        # Calculate insertion maneuver
        insertion_data = self.calculate_hyperbolic_insertion_maneuver(hyperbolic_orbit_data)
        
        # Create deployment schedule for all carriers
        deployment_schedules = self.create_deployment_schedule()
        
        # Display complete analysis
        results = self.display_mission_analysis(insertion_data, deployment_schedules)
        
        return results


def run_from_transfer_simulation():
    """
    Run deployment analysis using data from Transfer_Simulation
    """
    try:
        # Import Transfer_Simulation modules
        from core.interface import get_user_input
        from operations.earth_operations import calculate_earth_departure_delta_v
        from operations.trajectory_calculations import lunar_trajectory_calculations
        from operations.lunar_operations import lunar_soi_calculations, hyperbolic_to_elliptical_conversion
        
        print("Getting hyperbolic approach orbit data from Transfer_Simulation...")
        
        # Get transfer trajectory data
        user_inputs = get_user_input()
        R0 = user_inputs['R0']
        V0 = user_inputs['V0']
        gamma0 = user_inputs['gamma0']
        lambda1 = user_inputs['lambda1']
        
        # Calculate trajectory
        geo_results = lunar_trajectory_calculations(R0, V0, gamma0, lambda1, verbose=False)
        lunar_results = lunar_soi_calculations(
            geo_results['r1'], 
            geo_results['v1'], 
            geo_results['phi1_deg'],
            lambda1,
            geo_results['gamma1_deg'],
            verbose=False
        )
        
        # Extract hyperbolic data
        hyperbolic_data = {
            'v_hyp_perigee': math.sqrt((1 + lunar_results['e_lunar']) * MU_MOON_KMS / lunar_results['rp']),
            'perigee_altitude': lunar_results['rp'] - MOON_RADIUS_KM,
            'eccentricity': lunar_results['e_lunar']
        }
        
        print(f"Hyperbolic orbit data retrieved:")
        print(f"  Velocity at perigee: {hyperbolic_data['v_hyp_perigee']:.3f} km/s")
        print(f"  Perigee altitude: {hyperbolic_data['perigee_altitude']:.0f} km")
        print(f"  Eccentricity: {hyperbolic_data['eccentricity']:.3f}")
        
        # Run deployment analysis
        deployer = MoonConstellationDeployer()
        results = deployer.run_deployment_analysis(hyperbolic_data)
        
        return results
        
    except ImportError:
        print("Transfer_Simulation modules not available. Using default hyperbolic data.")
        deployer = MoonConstellationDeployer()
        results = deployer.run_deployment_analysis()
        return results


def main():
    """Main function"""
    print("Select deployment analysis mode:")
    print("1. Use Transfer_Simulation data (interactive)")
    print("2. Use default hyperbolic orbit data")
    
    try:
        choice = input("Enter choice (1 or 2): ").strip()
        
        if choice == '1':
            results = run_from_transfer_simulation()
        else:
            deployer = MoonConstellationDeployer()
            results = deployer.run_deployment_analysis()
            
        return results
        
    except KeyboardInterrupt:
        print("\nOperation cancelled by user")
        return None
    except Exception as e:
        print(f"Error: {e}")
        # Fallback to default analysis
        deployer = MoonConstellationDeployer()
        results = deployer.run_deployment_analysis()
        return results


if __name__ == "__main__":
    results = main()
