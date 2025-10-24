"""
Satellite Deployment System for Moon Communications Constellation
Starting from Elliptical 123km x 1500km Orbit at 62.8° Inclination

This system implements the deployment strategy:
1. Start from elliptical 123km x 1500km orbit at 62.8° inclination
2. Circularize at apoapsis to create 1500km circular orbit
3. Deploy 7 carriers (6 satellites each) with equal time intervals
4. Each carrier performs inclination change to 90° polar orbit
5. Each carrier deploys 6 satellites staggered to 1000km orbit

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
R_PERIGEE = MOON_RADIUS_KM + 123    # Initial perigee (123 km altitude)
R_APOGEE = MOON_RADIUS_KM + 1500    # Initial apogee (1500 km altitude)
R_CIRCULAR = MOON_RADIUS_KM + 1500  # Circularized orbit (1500 km altitude)
R_FINAL = MOON_RADIUS_KM + 1000     # Final satellite orbit (1000 km altitude)
INITIAL_INCLINATION = 62.8          # Initial inclination (degrees)
TARGET_INCLINATION = 90.0           # Target polar inclination (degrees)
NUM_CARRIERS = 7                    # Number of carriers
SATELLITES_PER_CARRIER = 6          # Satellites per carrier
TARGET_SPACING_DEG = 60             # degrees between satellites within each carrier
TOTAL_SATELLITES = NUM_CARRIERS * SATELLITES_PER_CARRIER  # 42 total satellites

class MoonConstellationDeployer:
    """
    Handles deployment of 42 satellites (7 carriers × 6 satellites each) starting from elliptical orbit:
    1. Start from elliptical 123km × 1500km orbit at 62.8° inclination
    2. Circularize at apoapsis to 1500km circular orbit
    3. Deploy 7 carriers with equal time intervals in 1500km orbit
    4. Each carrier performs inclination change to 90° polar orbit
    5. Each carrier deploys 6 satellites staggered to 1000km orbit
    
    Supports two deployment scenarios:
    - Single launcher: One spacecraft deploys all 7 carriers
    - Multi-launcher: 3 launchers (2+2+3 carriers) with coordinated deployment
    """
    
    def __init__(self, multi_launcher=False):
        self.moon_radius = MOON_RADIUS_KM
        self.mu_moon = MU_MOON_KMS
        
        # Orbital radii
        self.r_perigee = R_PERIGEE      # Initial perigee radius (123 km altitude)
        self.r_apogee = R_APOGEE        # Initial apogee radius (1500 km altitude)
        self.r_circular = R_CIRCULAR    # Circularized orbit radius (1500 km altitude)
        self.r_final = R_FINAL          # Final satellite orbit radius (1000 km altitude)
        
        # Mission parameters
        self.initial_inclination = INITIAL_INCLINATION
        self.target_inclination = TARGET_INCLINATION
        self.num_carriers = NUM_CARRIERS
        self.satellites_per_carrier = SATELLITES_PER_CARRIER
        self.target_spacing_deg = TARGET_SPACING_DEG
        self.total_satellites = TOTAL_SATELLITES
        
        # Multi-launcher configuration
        self.multi_launcher = multi_launcher
        if multi_launcher:
            self.launcher_config = [
                {'launcher_id': 1, 'num_carriers': 2, 'carrier_slots': [1, 2]},
                {'launcher_id': 2, 'num_carriers': 2, 'carrier_slots': [3, 4]},
                {'launcher_id': 3, 'num_carriers': 3, 'carrier_slots': [5, 6, 7]}
            ]
            self.num_launchers = len(self.launcher_config)
        else:
            self.launcher_config = [{'launcher_id': 1, 'num_carriers': 7, 'carrier_slots': list(range(1, 8))}]
            self.num_launchers = 1
        
        # Calculate all orbital parameters
        self.calculate_initial_elliptical_orbit()
        self.calculate_circularization_parameters()
        self.calculate_inclination_change_parameters()
        self.calculate_satellite_transfer_parameters()
        self.calculate_satellite_staggering_time()
    
    def calculate_initial_elliptical_orbit(self):
        """Calculate parameters for initial elliptical orbit"""
        # Semi-major axis and eccentricity
        self.a_initial = (self.r_perigee + self.r_apogee) / 2
        self.e_initial = (self.r_apogee - self.r_perigee) / (self.r_apogee + self.r_perigee)
        
        # Velocities at perigee and apogee
        self.v_perigee = math.sqrt(self.mu_moon * (2/self.r_perigee - 1/self.a_initial))
        self.v_apogee = math.sqrt(self.mu_moon * (2/self.r_apogee - 1/self.a_initial))
        
        # Orbital period
        self.T_initial = 2 * math.pi * math.sqrt(self.a_initial**3 / self.mu_moon)
        self.T_initial_hours = self.T_initial / 3600
        
        print(f"Initial Elliptical Orbit Parameters:")
        print(f"  Perigee: {self.r_perigee - self.moon_radius:.1f} km altitude")
        print(f"  Apogee: {self.r_apogee - self.moon_radius:.1f} km altitude")
        print(f"  Semi-major axis: {self.a_initial:.1f} km")
        print(f"  Eccentricity: {self.e_initial:.4f}")
        print(f"  Inclination: {self.initial_inclination:.1f}°")
        print(f"  Velocity at perigee: {self.v_perigee:.3f} km/s")
        print(f"  Velocity at apogee: {self.v_apogee:.3f} km/s")
        print(f"  Orbital period: {self.T_initial_hours:.2f} hours")
    
    def calculate_circularization_parameters(self):
        """Calculate parameters for circularization at apogee"""
        # Circular velocity at 1500km
        self.v_circular_1500 = math.sqrt(self.mu_moon / self.r_circular)
        
        # Delta-V for circularization (performed at apogee)
        self.delta_v_circularization = self.v_circular_1500 - self.v_apogee
        
        # Circular orbit period at 1500km
        self.T_circular = 2 * math.pi * math.sqrt(self.r_circular**3 / self.mu_moon)
        self.T_circular_hours = self.T_circular / 3600
        
        print(f"\nCircularization at Apogee (1500km):")
        print(f"  Required circular velocity: {self.v_circular_1500:.3f} km/s")
        print(f"  Current velocity at apogee: {self.v_apogee:.3f} km/s")
        print(f"  Delta-V for circularization: {self.delta_v_circularization*1000:.2f} m/s")
        print(f"  Circular orbit period: {self.T_circular_hours:.2f} hours")
    
    def calculate_inclination_change_parameters(self):
        """Calculate parameters for inclination change from 62.8° to 90°"""
        # Inclination change
        self.inclination_change = self.target_inclination - self.initial_inclination
        self.inclination_change_rad = math.radians(self.inclination_change)
        
        # Delta-V for inclination change (performed in circular 1500km orbit)
        # Formula: ΔV = 2 * V * sin(Δi/2)
        self.delta_v_inclination = 2 * self.v_circular_1500 * math.sin(self.inclination_change_rad / 2)
        
        print(f"\nInclination Change (per carrier):")
        print(f"  Initial inclination: {self.initial_inclination:.1f}°")
        print(f"  Target inclination: {self.target_inclination:.1f}°")
        print(f"  Inclination change: {self.inclination_change:.1f}°")
        print(f"  Orbital velocity (1500km): {self.v_circular_1500:.3f} km/s")
        print(f"  Delta-V for inclination change: {self.delta_v_inclination*1000:.2f} m/s")
    
    def calculate_satellite_transfer_parameters(self):
        """Calculate Hohmann transfer parameters from 1500km to 1000km"""
        # Transfer ellipse parameters
        self.a_sat_transfer = (self.r_circular + self.r_final) / 2  # Semi-major axis
        
        # Velocities for satellite Hohmann transfer
        # At 1500 km (apogee of transfer orbit)
        self.v_sat_transfer_1500 = math.sqrt(self.mu_moon * (2/self.r_circular - 1/self.a_sat_transfer))
        
        # At 1000 km (perigee of transfer orbit)  
        self.v_sat_transfer_1000 = math.sqrt(self.mu_moon * (2/self.r_final - 1/self.a_sat_transfer))
        
        # Circular velocities
        self.v_final_circular = math.sqrt(self.mu_moon / self.r_final)
        
        # Delta-V calculations for satellite transfer
        self.delta_v_sat1 = abs(self.v_circular_1500 - self.v_sat_transfer_1500)  # First burn (retrograde at 1500 km)
        self.delta_v_sat2 = abs(self.v_final_circular - self.v_sat_transfer_1000)  # Second burn (prograde at 1000 km)
        self.total_delta_v_sat = self.delta_v_sat1 + self.delta_v_sat2
        
        # Transfer time (half period of transfer ellipse)
        self.t_sat_transfer = math.pi * math.sqrt(self.a_sat_transfer**3 / self.mu_moon)
        self.t_sat_transfer_hours = self.t_sat_transfer / 3600
        
        # Final orbit parameters
        self.T_final = 2 * math.pi * math.sqrt(self.r_final**3 / self.mu_moon)
        self.T_final_hours = self.T_final / 3600
        self.n_final = 2 * math.pi / self.T_final  # Mean motion final orbit
        self.n_circular = 2 * math.pi / self.T_circular  # Mean motion 1500km orbit
        
        print(f"\nSatellite Hohmann Transfer (1500km → 1000km):")
        print(f"  Transfer semi-major axis: {self.a_sat_transfer:.1f} km")
        print(f"  First burn (retrograde at 1500km): {self.delta_v_sat1*1000:.2f} m/s")
        print(f"  Second burn (prograde at 1000km): {self.delta_v_sat2*1000:.2f} m/s")
        print(f"  Total transfer Delta-V per satellite: {self.total_delta_v_sat*1000:.2f} m/s")
        print(f"  Transfer time: {self.t_sat_transfer_hours:.3f} hours")
        print(f"  Final orbit period (1000km): {self.T_final_hours:.2f} hours")
    
    def calculate_satellite_staggering_time(self):
        """Calculate the staggering time τ for 60° spacing in final 1000km orbit"""
        # The staggering time is based on differential motion between 1500km and 1000km orbits
        # Δθₖ = (n_final - n_circular)(k-1)τ
        # For 60° spacing: Δθ = π/3 radians
        # τ = π/3 / (n_final - n_circular)
        
        spacing_rad = math.radians(self.target_spacing_deg)
        self.tau_sat = spacing_rad / (self.n_final - self.n_circular)
        self.tau_sat_hours = self.tau_sat / 3600
        
        print(f"\nSatellite Staggering Parameters:")
        print(f"  Target spacing: {self.target_spacing_deg}° = {spacing_rad:.4f} rad")
        print(f"  Mean motion difference: n_final - n_circular = {(self.n_final - self.n_circular):.4e} rad/s")
        print(f"  Staggering time τ = {self.tau_sat:.1f} s = {self.tau_sat_hours:.3f} hours")
    
    def calculate_carrier_deployment_schedule(self):
        """Calculate deployment schedule for 7 carriers"""
        # Time interval between carrier deployments (equal spacing in 1500km orbit)
        carrier_spacing_time = self.T_circular / self.num_carriers
        carrier_spacing_hours = carrier_spacing_time / 3600
        
        print(f"\nCarrier Deployment Schedule:")
        print(f"  Number of carriers: {self.num_carriers}")
        print(f"  Deployment orbit period: {self.T_circular_hours:.2f} hours")
        print(f"  Time interval between carriers: {carrier_spacing_hours:.3f} hours ({carrier_spacing_time:.1f} s)")
        
        return carrier_spacing_time, carrier_spacing_hours
    
    def calculate_multi_launcher_delta_v_budget(self):
        """Calculate total Delta-V budget for multi-launcher mission"""
        
        # Each launcher performs circularization
        launcher_circularization = self.delta_v_circularization * 1000  # m/s per launcher
        total_launcher_delta_v = self.num_launchers * launcher_circularization
        
        # Carrier Delta-V (per carrier) - same as single launcher
        carrier_inclination_change = self.delta_v_inclination * 1000  # m/s per carrier
        total_carrier_delta_v = self.num_carriers * carrier_inclination_change
        
        # Satellite Delta-V (per satellite) - same as single launcher
        satellite_transfer = self.total_delta_v_sat * 1000  # m/s per satellite
        total_satellite_delta_v = self.total_satellites * satellite_transfer
        
        # Total mission Delta-V
        total_mission_delta_v = total_launcher_delta_v + total_carrier_delta_v + total_satellite_delta_v
        
        delta_v_budget = {
            'launchers': {
                'circularization_per_launcher_ms': launcher_circularization,
                'num_launchers': self.num_launchers,
                'total_ms': total_launcher_delta_v
            },
            'carriers': {
                'inclination_change_per_carrier_ms': carrier_inclination_change,
                'num_carriers': self.num_carriers,
                'total_ms': total_carrier_delta_v
            },
            'satellites': {
                'transfer_per_satellite_ms': satellite_transfer,
                'num_satellites': self.total_satellites,
                'total_ms': total_satellite_delta_v
            },
            'mission_total': {
                'total_ms': total_mission_delta_v,
                'total_kms': total_mission_delta_v / 1000
            }
        }
        
        return delta_v_budget
    
    def create_multi_launcher_mission_schedule(self):
        """Create complete mission schedule for multi-launcher scenario"""
        
        # All launchers start simultaneously and circularize
        circularization_time = 0  # All performed simultaneously at first apogee passage
        
        # Calculate carrier deployment timing based on predetermined slots
        carrier_spacing_time, carrier_spacing_hours = self.calculate_carrier_deployment_schedule()
        
        mission_schedule = {
            'phase1_circularization': {
                'description': 'All 3 launchers circularize orbits at 1500km apogee',
                'launchers': []
            },
            'phase2_carrier_deployment': {
                'description': 'Coordinated carrier deployment following slot assignments',
                'carriers': []
            },
            'phase3_inclination_changes': {
                'description': 'Each carrier changes inclination to 90°',
                'carriers': []
            },
            'phase4_satellite_deployment': {
                'description': 'Each carrier deploys 6 satellites to 1000km orbit',
                'carriers': []
            }
        }
        
        # Phase 1: All launchers circularize simultaneously
        for launcher_config in self.launcher_config:
            launcher_data = {
                'launcher_id': launcher_config['launcher_id'],
                'num_carriers': launcher_config['num_carriers'],
                'circularization_time_s': circularization_time,
                'circularization_time_h': circularization_time / 3600,
                'delta_v_ms': self.delta_v_circularization * 1000,
                'result_orbit': '1500km circular at 62.8° inclination'
            }
            mission_schedule['phase1_circularization']['launchers'].append(launcher_data)
        
        # Calculate carrier deployment schedule following slot assignments
        stabilization_time = self.T_circular  # One orbit for stabilization after circularization
        
        for launcher_config in self.launcher_config:
            launcher_id = launcher_config['launcher_id']
            
            for i, carrier_slot in enumerate(launcher_config['carrier_slots']):
                # Carrier deployment time based on predetermined slot
                carrier_deployment_time = stabilization_time + (carrier_slot - 1) * carrier_spacing_time
                
                # Inclination change performed immediately after deployment
                inclination_change_time = carrier_deployment_time + 300  # 5 minutes after deployment
                
                # Carrier data for phase 2
                carrier_deployment = {
                    'carrier_id': carrier_slot,
                    'launcher_id': launcher_id,
                    'deployment_time_s': carrier_deployment_time,
                    'deployment_time_h': carrier_deployment_time / 3600,
                    'deployment_slot': carrier_slot,
                    'deployment_orbit': '1500km circular at 62.8° inclination'
                }
                mission_schedule['phase2_carrier_deployment']['carriers'].append(carrier_deployment)
                
                # Carrier data for phase 3
                carrier_inclination = {
                    'carrier_id': carrier_slot,
                    'launcher_id': launcher_id,
                    'inclination_change_time_s': inclination_change_time,
                    'inclination_change_time_h': inclination_change_time / 3600,
                    'delta_v_ms': self.delta_v_inclination * 1000,
                    'initial_inclination_deg': self.initial_inclination,
                    'final_inclination_deg': self.target_inclination,
                    'result_orbit': '1500km circular at 90° inclination (polar)'
                }
                mission_schedule['phase3_inclination_changes']['carriers'].append(carrier_inclination)
                
                # Satellite deployment for this carrier (phase 4)
                satellite_deployment_start = inclination_change_time + self.T_circular  # One orbit after inclination change
                
                carrier_satellites = {
                    'carrier_id': carrier_slot,
                    'launcher_id': launcher_id,
                    'satellites': []
                }
                
                # Calculate staggered satellite deployments for this carrier
                for sat_id in range(1, self.satellites_per_carrier + 1):
                    sat_transfer_start = satellite_deployment_start + (sat_id - 1) * self.tau_sat
                    sat_first_burn = sat_transfer_start
                    sat_second_burn = sat_transfer_start + self.t_sat_transfer
                    
                    # Final angular position due to staggered transfers
                    final_angular_position_deg = (sat_id - 1) * self.target_spacing_deg
                    
                    satellite_data = {
                        'satellite_id': sat_id,
                        'global_satellite_id': (carrier_slot - 1) * self.satellites_per_carrier + sat_id,
                        'transfer_start_time_s': sat_transfer_start,
                        'transfer_start_time_h': sat_transfer_start / 3600,
                        'first_burn_time_s': sat_first_burn,
                        'first_burn_time_h': sat_first_burn / 3600,
                        'second_burn_time_s': sat_second_burn,
                        'second_burn_time_h': sat_second_burn / 3600,
                        'delta_v1_ms': self.delta_v_sat1 * 1000,
                        'delta_v2_ms': self.delta_v_sat2 * 1000,
                        'total_transfer_delta_v_ms': self.total_delta_v_sat * 1000,
                        'final_altitude_km': 1000,
                        'final_inclination_deg': 90,
                        'final_angular_position_deg': final_angular_position_deg
                    }
                    
                    carrier_satellites['satellites'].append(satellite_data)
                
                mission_schedule['phase4_satellite_deployment']['carriers'].append(carrier_satellites)
        
        return mission_schedule
    
    def create_complete_mission_schedule(self):
        """Create complete mission schedule from elliptical orbit to final constellation"""
        
        # Phase 1: Circularization at apogee
        circularization_time = 0  # Performed immediately at first apogee passage
        
        # Phase 2: Carrier deployment
        carrier_spacing_time, carrier_spacing_hours = self.calculate_carrier_deployment_schedule()
        
        mission_schedule = {
            'phase1_circularization': {
                'description': 'Circularize orbit at 1500km apogee',
                'time_s': circularization_time,
                'time_h': circularization_time / 3600,
                'delta_v_ms': self.delta_v_circularization * 1000,
                'result_orbit': '1500km circular at 62.8° inclination'
            },
            'phase2_carrier_deployment': {
                'description': 'Deploy 7 carriers with equal time intervals',
                'carriers': []
            },
            'phase3_inclination_changes': {
                'description': 'Each carrier changes inclination to 90°',
                'carriers': []
            },
            'phase4_satellite_deployment': {
                'description': 'Each carrier deploys 6 satellites to 1000km orbit',
                'carriers': []
            }
        }
        
        # Calculate deployment schedule for each carrier
        stabilization_time = self.T_circular  # One orbit for stabilization after circularization
        
        for carrier_id in range(1, self.num_carriers + 1):
            # Carrier deployment time
            carrier_deployment_time = stabilization_time + (carrier_id - 1) * carrier_spacing_time
            
            # Inclination change time (performed immediately after deployment)
            inclination_change_time = carrier_deployment_time + 300  # 5 minutes after deployment
            
            # Carrier data for phase 2
            carrier_deployment = {
                'carrier_id': carrier_id,
                'deployment_time_s': carrier_deployment_time,
                'deployment_time_h': carrier_deployment_time / 3600,
                'deployment_orbit': '1500km circular at 62.8° inclination'
            }
            mission_schedule['phase2_carrier_deployment']['carriers'].append(carrier_deployment)
            
            # Carrier data for phase 3
            carrier_inclination = {
                'carrier_id': carrier_id,
                'inclination_change_time_s': inclination_change_time,
                'inclination_change_time_h': inclination_change_time / 3600,
                'delta_v_ms': self.delta_v_inclination * 1000,
                'initial_inclination_deg': self.initial_inclination,
                'final_inclination_deg': self.target_inclination,
                'result_orbit': '1500km circular at 90° inclination (polar)'
            }
            mission_schedule['phase3_inclination_changes']['carriers'].append(carrier_inclination)
            
            # Satellite deployment for this carrier (phase 4)
            satellite_deployment_start = inclination_change_time + self.T_circular  # One orbit after inclination change
            
            carrier_satellites = {
                'carrier_id': carrier_id,
                'satellites': []
            }
            
            # Calculate staggered satellite deployments for this carrier
            for sat_id in range(1, self.satellites_per_carrier + 1):
                sat_transfer_start = satellite_deployment_start + (sat_id - 1) * self.tau_sat
                sat_first_burn = sat_transfer_start
                sat_second_burn = sat_transfer_start + self.t_sat_transfer
                
                # Final angular position due to staggered transfers
                final_angular_position_deg = (sat_id - 1) * self.target_spacing_deg
                
                satellite_data = {
                    'satellite_id': sat_id,
                    'global_satellite_id': (carrier_id - 1) * self.satellites_per_carrier + sat_id,
                    'transfer_start_time_s': sat_transfer_start,
                    'transfer_start_time_h': sat_transfer_start / 3600,
                    'first_burn_time_s': sat_first_burn,
                    'first_burn_time_h': sat_first_burn / 3600,
                    'second_burn_time_s': sat_second_burn,
                    'second_burn_time_h': sat_second_burn / 3600,
                    'delta_v1_ms': self.delta_v_sat1 * 1000,
                    'delta_v2_ms': self.delta_v_sat2 * 1000,
                    'total_transfer_delta_v_ms': self.total_delta_v_sat * 1000,
                    'final_altitude_km': 1000,
                    'final_inclination_deg': 90,
                    'final_angular_position_deg': final_angular_position_deg
                }
                
                carrier_satellites['satellites'].append(satellite_data)
            
            mission_schedule['phase4_satellite_deployment']['carriers'].append(carrier_satellites)
        
        return mission_schedule
    
    def calculate_total_delta_v_budget(self):
        """Calculate total Delta-V budget for the complete mission"""
        
        # Spacecraft Delta-V (main spacecraft)
        spacecraft_circularization = self.delta_v_circularization * 1000  # m/s
        
        # Carrier Delta-V (per carrier)
        carrier_inclination_change = self.delta_v_inclination * 1000  # m/s per carrier
        total_carrier_delta_v = self.num_carriers * carrier_inclination_change
        
        # Satellite Delta-V (per satellite)
        satellite_transfer = self.total_delta_v_sat * 1000  # m/s per satellite
        total_satellite_delta_v = self.total_satellites * satellite_transfer
        
        # Total mission Delta-V
        total_mission_delta_v = spacecraft_circularization + total_carrier_delta_v + total_satellite_delta_v
        
        delta_v_budget = {
            'spacecraft': {
                'circularization_ms': spacecraft_circularization,
                'total_ms': spacecraft_circularization
            },
            'carriers': {
                'inclination_change_per_carrier_ms': carrier_inclination_change,
                'num_carriers': self.num_carriers,
                'total_ms': total_carrier_delta_v
            },
            'satellites': {
                'transfer_per_satellite_ms': satellite_transfer,
                'num_satellites': self.total_satellites,
                'total_ms': total_satellite_delta_v
            },
            'mission_total': {
                'total_ms': total_mission_delta_v,
                'total_kms': total_mission_delta_v / 1000
            }
        }
        
        return delta_v_budget
    
    def display_multi_launcher_mission_analysis(self):
        """Display comprehensive mission analysis for multi-launcher scenario"""
        print("\n" + "=" * 100)
        print("MOON COMMUNICATIONS CONSTELLATION DEPLOYMENT ANALYSIS")
        print("MULTI-LAUNCHER SCENARIO: 3 Launchers with Coordinated Deployment")
        print("Starting from Elliptical 123km × 1500km Orbit at 62.8° Inclination")
        print("=" * 100)
        
        # Create mission schedule and calculate Delta-V budget
        mission_schedule = self.create_multi_launcher_mission_schedule()
        delta_v_budget = self.calculate_multi_launcher_delta_v_budget()
        
        print(f"\nMISSION OVERVIEW:")
        print(f"   Total satellites: {self.total_satellites} ({self.num_carriers} carriers × {self.satellites_per_carrier} satellites each)")
        print(f"   Number of launchers: {self.num_launchers}")
        print(f"   Launcher configuration: 2+2+3 carriers")
        print(f"   Initial orbit: 123km × 1500km elliptical at {self.initial_inclination}° inclination")
        print(f"   Intermediate orbit: 1500km circular at {self.initial_inclination}° inclination")
        print(f"   Carrier final orbits: 1500km circular at {self.target_inclination}° inclination (polar)")
        print(f"   Satellite final orbits: 1000km circular at {self.target_inclination}° inclination (polar)")
        
        # Display launcher configuration
        print(f"\nLAUNCHER CONFIGURATION:")
        for launcher_config in self.launcher_config:
            launcher_id = launcher_config['launcher_id']
            num_carriers = launcher_config['num_carriers']
            carrier_slots = launcher_config['carrier_slots']
            print(f"   Launcher {launcher_id}: {num_carriers} carriers (slots {carrier_slots})")
        
        # Phase 1: Circularization
        print(f"\nPHASE 1: ORBIT CIRCULARIZATION (ALL LAUNCHERS)")
        print(f"   Initial orbit: {self.r_perigee - self.moon_radius:.0f}km × {self.r_apogee - self.moon_radius:.0f}km elliptical")
        print(f"   Maneuver: Prograde burn at apogee (1500km altitude)")
        print(f"   Delta-V per launcher: {delta_v_budget['launchers']['circularization_per_launcher_ms']:.1f} m/s")
        print(f"   Total for all launchers: {delta_v_budget['launchers']['total_ms']:.1f} m/s")
        print(f"   Result: 1500km circular orbit at {self.initial_inclination}° inclination")
        print(f"   Time: All performed simultaneously at t = 0")
        
        # Phase 2: Carrier deployment
        print(f"\nPHASE 2: COORDINATED CARRIER DEPLOYMENT")
        print(f"   Total carriers: {self.num_carriers}")
        print(f"   Deployment orbit: 1500km circular at {self.initial_inclination}° inclination")
        print(f"   Deployment method: Coordinated timing following predetermined slots")
        carrier_spacing_time, carrier_spacing_hours = self.calculate_carrier_deployment_schedule()
        print(f"   Time interval between carrier slots: {carrier_spacing_hours:.3f} hours ({carrier_spacing_time:.1f} seconds)")
        
        print(f"\n   Carrier Deployment Timeline:")
        print(f"   {'Launcher':<10} {'Carrier':<8} {'Slot':<6} {'Time (h)':<10} {'Carriers':<15}")
        print(f"   {'-'*10} {'-'*8} {'-'*6} {'-'*10} {'-'*15}")
        
        for carrier in sorted(mission_schedule['phase2_carrier_deployment']['carriers'], key=lambda x: x['deployment_slot']):
            launcher_id = carrier['launcher_id']
            carrier_id = carrier['carrier_id']
            slot = carrier['deployment_slot']
            time_h = carrier['deployment_time_h']
            
            # Find launcher config for this launcher
            launcher_config = next(lc for lc in self.launcher_config if lc['launcher_id'] == launcher_id)
            carrier_list = f"[{', '.join(map(str, launcher_config['carrier_slots']))}]"
            
            print(f"   {launcher_id:<10} {carrier_id:<8} {slot:<6} {time_h:<10.3f} {carrier_list:<15}")
        
        print(f"   Total deployment span: {(self.num_carriers - 1) * carrier_spacing_hours:.2f} hours")
        
        # Phase 3: Inclination changes
        print(f"\nPHASE 3: INCLINATION CHANGES (PER CARRIER)")
        print(f"   Maneuver: Change inclination from {self.initial_inclination}° to {self.target_inclination}°")
        print(f"   Delta-V per carrier: {delta_v_budget['carriers']['inclination_change_per_carrier_ms']:.1f} m/s")
        print(f"   Total for all carriers: {delta_v_budget['carriers']['total_ms']:.1f} m/s")
        print(f"   Orbit after maneuver: 1500km circular polar orbit")
        
        # Phase 4: Satellite deployment
        print(f"\nPHASE 4: SATELLITE DEPLOYMENT (PER CARRIER)")
        print(f"   Satellites per carrier: {self.satellites_per_carrier}")
        print(f"   Transfer method: Staggered Hohmann transfers (1500km → 1000km)")
        print(f"   Staggering time: {self.tau_sat_hours:.3f} hours ({self.tau_sat:.1f} seconds)")
        print(f"   Transfer time per satellite: {self.t_sat_transfer_hours:.3f} hours")
        print(f"   Delta-V per satellite: {delta_v_budget['satellites']['transfer_per_satellite_ms']:.1f} m/s")
        print(f"   Total for all satellites: {delta_v_budget['satellites']['total_ms']:.1f} m/s")
        print(f"   Final spacing: {self.target_spacing_deg}° between satellites within each carrier")
        
        # Calculate total mission timeline
        last_satellite_completion = 0
        for carrier_data in mission_schedule['phase4_satellite_deployment']['carriers']:
            for sat_data in carrier_data['satellites']:
                completion_time = sat_data['second_burn_time_h']
                if completion_time > last_satellite_completion:
                    last_satellite_completion = completion_time
        
        print(f"\nDELTA-V BUDGET SUMMARY (MULTI-LAUNCHER):")
        print(f"   {'Component':<25} {'Delta-V (m/s)':<15} {'Count':<10} {'Total (m/s)':<15}")
        print(f"   {'-'*25} {'-'*15} {'-'*10} {'-'*15}")
        print(f"   {'Launcher circularization':<25} {delta_v_budget['launchers']['circularization_per_launcher_ms']:<15.1f} {delta_v_budget['launchers']['num_launchers']:<10} {delta_v_budget['launchers']['total_ms']:<15.1f}")
        print(f"   {'Carrier inclination change':<25} {delta_v_budget['carriers']['inclination_change_per_carrier_ms']:<15.1f} {delta_v_budget['carriers']['num_carriers']:<10} {delta_v_budget['carriers']['total_ms']:<15.1f}")
        print(f"   {'Satellite transfers':<25} {delta_v_budget['satellites']['transfer_per_satellite_ms']:<15.1f} {delta_v_budget['satellites']['num_satellites']:<10} {delta_v_budget['satellites']['total_ms']:<15.1f}")
        print(f"   {'-'*25} {'-'*15} {'-'*10} {'-'*15}")
        print(f"   {'TOTAL MISSION':<25} {'':<15} {'':<10} {delta_v_budget['mission_total']['total_ms']:<15.1f}")
        print(f"   {'TOTAL MISSION (km/s)':<25} {'':<15} {'':<10} {delta_v_budget['mission_total']['total_kms']:<15.3f}")
        
        print(f"\nMISSION TIMELINE SUMMARY (MULTI-LAUNCHER):")
        print(f"   Phase 1 (Circularization): 0.0 hours (all launchers simultaneous)")
        print(f"   Phase 2 (Carrier deployment): {(self.num_carriers - 1) * carrier_spacing_hours:.2f} hours")
        print(f"   Phase 3 (Inclination changes): Overlapped with Phase 2")
        print(f"   Phase 4 (Satellite deployment): {last_satellite_completion:.2f} hours")
        print(f"   Total mission duration: {last_satellite_completion:.2f} hours ({last_satellite_completion/24:.1f} days)")
        
        print(f"\nMULTI-LAUNCHER ADVANTAGES:")
        print(f"   • Reduced risk through mission redundancy")
        print(f"   • Parallel operations increase deployment efficiency")
        print(f"   • Distributed payload mass across multiple launches")
        print(f"   • Coordinated timing maintains constellation geometry")
        print(f"   • Same total Delta-V as single launcher scenario")
        
        print(f"\nFINAL CONSTELLATION:")
        print(f"   Total satellites: {self.total_satellites}")
        print(f"   Orbital planes: {self.num_carriers} (polar orbits at different RAAN)")
        print(f"   Satellites per plane: {self.satellites_per_carrier}")
        print(f"   Satellite altitude: 1000km")
        print(f"   Satellite inclination: 90° (polar)")
        print(f"   Spacing within each plane: {self.target_spacing_deg}°")
        print(f"   Coverage: Global lunar coverage")
        print(f"   Constellation revisit time: ~{self.T_final_hours/self.satellites_per_carrier:.2f} hours per plane")
        
        return {
            'mission_schedule': mission_schedule,
            'delta_v_budget': delta_v_budget,
            'total_mission_time_h': last_satellite_completion,
            'total_mission_delta_v_ms': delta_v_budget['mission_total']['total_ms']
        }

    def display_mission_analysis(self):
        """Display comprehensive mission analysis"""
        if self.multi_launcher:
            return self.display_multi_launcher_mission_analysis()
        
        print("\n" + "=" * 100)
        print("MOON COMMUNICATIONS CONSTELLATION DEPLOYMENT ANALYSIS")
        print("SINGLE LAUNCHER SCENARIO")
        print("Starting from Elliptical 123km × 1500km Orbit at 62.8° Inclination")
        print("=" * 100)
        
        # Create mission schedule and calculate Delta-V budget
        mission_schedule = self.create_complete_mission_schedule()
        delta_v_budget = self.calculate_total_delta_v_budget()
        
        print(f"\nMISSION OVERVIEW:")
        print(f"   Total satellites: {self.total_satellites} ({self.num_carriers} carriers × {self.satellites_per_carrier} satellites each)")
        print(f"   Initial orbit: 123km × 1500km elliptical at {self.initial_inclination}° inclination")
        print(f"   Intermediate orbit: 1500km circular at {self.initial_inclination}° inclination")
        print(f"   Carrier final orbits: 1500km circular at {self.target_inclination}° inclination (polar)")
        print(f"   Satellite final orbits: 1000km circular at {self.target_inclination}° inclination (polar)")
        
        # Phase 1: Circularization
        print(f"\nPHASE 1: ORBIT CIRCULARIZATION")
        print(f"   Initial orbit: {self.r_perigee - self.moon_radius:.0f}km × {self.r_apogee - self.moon_radius:.0f}km elliptical")
        print(f"   Maneuver: Prograde burn at apogee (1500km altitude)")
        print(f"   Delta-V required: {delta_v_budget['spacecraft']['circularization_ms']:.1f} m/s")
        print(f"   Result: 1500km circular orbit at {self.initial_inclination}° inclination")
        print(f"   Time: Performed at first apogee passage (t = 0)")
        
        # Phase 2: Carrier deployment
        print(f"\nPHASE 2: CARRIER DEPLOYMENT")
        print(f"   Number of carriers: {self.num_carriers}")
        print(f"   Deployment orbit: 1500km circular at {self.initial_inclination}° inclination")
        print(f"   Deployment method: Equal time intervals")
        carrier_spacing_time, carrier_spacing_hours = self.calculate_carrier_deployment_schedule()
        print(f"   Time interval between carriers: {carrier_spacing_hours:.3f} hours ({carrier_spacing_time:.1f} seconds)")
        print(f"   Total deployment time: {(self.num_carriers - 1) * carrier_spacing_hours:.2f} hours")
        
        # Phase 3: Inclination changes
        print(f"\nPHASE 3: INCLINATION CHANGES (PER CARRIER)")
        print(f"   Maneuver: Change inclination from {self.initial_inclination}° to {self.target_inclination}°")
        print(f"   Delta-V per carrier: {delta_v_budget['carriers']['inclination_change_per_carrier_ms']:.1f} m/s")
        print(f"   Total for all carriers: {delta_v_budget['carriers']['total_ms']:.1f} m/s")
        print(f"   Orbit after maneuver: 1500km circular polar orbit")
        
        # Phase 4: Satellite deployment
        print(f"\nPHASE 4: SATELLITE DEPLOYMENT (PER CARRIER)")
        print(f"   Satellites per carrier: {self.satellites_per_carrier}")
        print(f"   Transfer method: Staggered Hohmann transfers (1500km → 1000km)")
        print(f"   Staggering time: {self.tau_sat_hours:.3f} hours ({self.tau_sat:.1f} seconds)")
        print(f"   Transfer time per satellite: {self.t_sat_transfer_hours:.3f} hours")
        print(f"   Delta-V per satellite: {delta_v_budget['satellites']['transfer_per_satellite_ms']:.1f} m/s")
        print(f"   Total for all satellites: {delta_v_budget['satellites']['total_ms']:.1f} m/s")
        print(f"   Final spacing: {self.target_spacing_deg}° between satellites within each carrier")
        
        # Calculate total mission timeline
        last_satellite_completion = 0
        for carrier_data in mission_schedule['phase4_satellite_deployment']['carriers']:
            for sat_data in carrier_data['satellites']:
                completion_time = sat_data['second_burn_time_h']
                if completion_time > last_satellite_completion:
                    last_satellite_completion = completion_time
        
        print(f"\nDELTA-V BUDGET SUMMARY:")
        print(f"   {'Component':<25} {'Delta-V (m/s)':<15} {'Count':<10} {'Total (m/s)':<15}")
        print(f"   {'-'*25} {'-'*15} {'-'*10} {'-'*15}")
        print(f"   {'Spacecraft circularization':<25} {delta_v_budget['spacecraft']['circularization_ms']:<15.1f} {'1':<10} {delta_v_budget['spacecraft']['total_ms']:<15.1f}")
        print(f"   {'Carrier inclination change':<25} {delta_v_budget['carriers']['inclination_change_per_carrier_ms']:<15.1f} {delta_v_budget['carriers']['num_carriers']:<10} {delta_v_budget['carriers']['total_ms']:<15.1f}")
        print(f"   {'Satellite transfers':<25} {delta_v_budget['satellites']['transfer_per_satellite_ms']:<15.1f} {delta_v_budget['satellites']['num_satellites']:<10} {delta_v_budget['satellites']['total_ms']:<15.1f}")
        print(f"   {'-'*25} {'-'*15} {'-'*10} {'-'*15}")
        print(f"   {'TOTAL MISSION':<25} {'':<15} {'':<10} {delta_v_budget['mission_total']['total_ms']:<15.1f}")
        print(f"   {'TOTAL MISSION (km/s)':<25} {'':<15} {'':<10} {delta_v_budget['mission_total']['total_kms']:<15.3f}")
        
        print(f"\nMISSION TIMELINE SUMMARY:")
        print(f"   Phase 1 (Circularization): 0.0 hours")
        print(f"   Phase 2 (Carrier deployment): {(self.num_carriers - 1) * carrier_spacing_hours:.2f} hours")
        print(f"   Phase 3 (Inclination changes): Overlapped with Phase 2")
        print(f"   Phase 4 (Satellite deployment): {last_satellite_completion:.2f} hours")
        print(f"   Total mission duration: {last_satellite_completion:.2f} hours ({last_satellite_completion/24:.1f} days)")
        
        print(f"\nFINAL CONSTELLATION:")
        print(f"   Total satellites: {self.total_satellites}")
        print(f"   Orbital planes: {self.num_carriers} (polar orbits at different RAAN)")
        print(f"   Satellites per plane: {self.satellites_per_carrier}")
        print(f"   Satellite altitude: 1000km")
        print(f"   Satellite inclination: 90° (polar)")
        print(f"   Spacing within each plane: {self.target_spacing_deg}°")
        print(f"   Coverage: Global lunar coverage")
        print(f"   Constellation revisit time: ~{self.T_final_hours/self.satellites_per_carrier:.2f} hours per plane")
        
        return {
            'mission_schedule': mission_schedule,
            'delta_v_budget': delta_v_budget,
            'total_mission_time_h': last_satellite_completion,
            'total_mission_delta_v_ms': delta_v_budget['mission_total']['total_ms']
        }
    
    def run_deployment_analysis(self):
        """Run complete constellation deployment analysis"""
        print("MOON CONSTELLATION DEPLOYER")
        print("Deployment from Elliptical 123km × 1500km Orbit at 62.8° Inclination")
        print("=" * 80)
        
        # Run complete analysis
        results = self.display_mission_analysis()
        
        return results


def run_single_launcher_analysis():
    """
    Run analysis for single launcher scenario
    """
    print("Starting single launcher constellation deployment...")
    
    deployer = MoonConstellationDeployer(multi_launcher=False)
    results = deployer.run_deployment_analysis()
    
    return results


def run_multi_launcher_analysis():
    """
    Run analysis for multi-launcher scenario (3 launchers: 2+2+3 carriers)
    """
    print("Starting multi-launcher constellation deployment...")
    
    deployer = MoonConstellationDeployer(multi_launcher=True)
    results = deployer.run_deployment_analysis()
    
    return results


def main():
    """Main function"""
    print("Moon Communications Constellation Deployment Analysis")
    print("Starting from Elliptical 123km × 1500km Orbit at 62.8° Inclination")
    print("=" * 80)
    print("\nSelect deployment scenario:")
    print("1. Single launcher (1 spacecraft deploys all 7 carriers)")
    print("2. Multi-launcher (3 launchers: 2+2+3 carriers)")
    print("3. Both scenarios (comprehensive comparison)")
    
    try:
        choice = input("\nEnter choice (1, 2, or 3): ").strip()
        
        if choice == '1':
            results = run_single_launcher_analysis()
        elif choice == '2':
            results = run_multi_launcher_analysis()
        elif choice == '3':
            print("\n" + "="*100)
            print("COMPREHENSIVE ANALYSIS: BOTH DEPLOYMENT SCENARIOS")
            print("="*100)
            
            print("\n" + "="*50)
            print("SCENARIO 1: SINGLE LAUNCHER")
            print("="*50)
            results_single = run_single_launcher_analysis()
            
            print("\n\n" + "="*50)  
            print("SCENARIO 2: MULTI-LAUNCHER")
            print("="*50)
            results_multi = run_multi_launcher_analysis()
            
            # Comparison summary
            print("\n" + "="*100)
            print("DEPLOYMENT SCENARIOS COMPARISON")
            print("="*100)
            
            print(f"\n{'Metric':<30} {'Single Launcher':<20} {'Multi-Launcher':<20}")
            print(f"{'-'*30} {'-'*20} {'-'*20}")
            print(f"{'Number of launchers':<30} {'1':<20} {'3':<20}")
            print(f"{'Total Delta-V (m/s)':<30} {results_single['total_mission_delta_v_ms']:<20.1f} {results_multi['total_mission_delta_v_ms']:<20.1f}")
            print(f"{'Mission duration (hours)':<30} {results_single['total_mission_time_h']:<20.2f} {results_multi['total_mission_time_h']:<20.2f}")
            print(f"{'Mission duration (days)':<30} {results_single['total_mission_time_h']/24:<20.1f} {results_multi['total_mission_time_h']/24:<20.1f}")
            
            results = {'single_launcher': results_single, 'multi_launcher': results_multi}
        else:
            print("Invalid choice. Running single launcher analysis...")
            results = run_single_launcher_analysis()
            
        return results
        
    except KeyboardInterrupt:
        print("\nOperation cancelled by user")
        return None
    except Exception as e:
        print(f"Error: {e}")
        return None


if __name__ == "__main__":
    results = main()
