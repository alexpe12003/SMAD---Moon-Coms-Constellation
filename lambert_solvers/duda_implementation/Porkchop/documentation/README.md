## Circular Orbit Position and Time Calculator - Summary

### Project Overview
This project successfully defines the position and time for two circular orbits:
1. **LEO Satellite**: 200km altitude circular orbit around Earth
2. **Moon**: Circular orbit at the average Earth-Moon distance (384,400 km)

Both orbits are coplanar and 2D, with position defined by angle θ where θ=0 at the initial position, calculated at 1-minute intervals for the duration of one complete LEO satellite orbit (88.5 minutes).

### Files Created

#### 1. `porkchop.py` - Main Orbital Calculator
- **CircularOrbit class**: Core orbital mechanics calculations
- **Visualization functions**: Plotting orbital trajectories and angular positions
- **Data generation**: Creates position arrays with customizable time intervals
- **Key Features**:
  - Calculates orbital velocity: v = √(μ/r)
  - Calculates orbital period: T = 2π√(r³/μ)
  - Generates position coordinates: x(t) = r·cos(θ(t)), y(t) = r·sin(θ(t))

#### 2. `detailed_analysis.py` - Mathematical Analysis
- **Mathematical derivations**: Step-by-step orbital mechanics calculations
- **Comparative analysis**: Ratio comparisons between LEO and Moon orbits
- **Synodic period calculations**: Time between successive alignments
- **Position equation demonstrations**: Examples at specific time points

#### 3. `data_access.py` - Data Access and Export
- **Data array generation**: Creates complete position/time datasets
- **Sample data display**: Shows first 20 data points in tabular format
- **Usage examples**: Demonstrates how to access and use the data arrays
- **File export**: Saves all data to `orbit_data.txt` for external use

#### 4. `orbit_data.txt` - Complete Dataset
- Contains 90 data points (88.5 minutes at 1-minute intervals)
- Columns: Time, LEO θ, LEO X, LEO Y, Moon θ, Moon X, Moon Y
- UTF-8 encoded text file ready for analysis in other applications

### Key Results

#### LEO Satellite (200km altitude)
- **Orbital radius**: 6,578 km
- **Orbital velocity**: 7.784 km/s
- **Orbital period**: 88.5 minutes (1.47 hours)
- **Angular velocity**: 0.068 deg/s

#### Moon Orbit
- **Orbital radius**: 384,400 km
- **Orbital velocity**: 1.018 km/s
- **Orbital period**: 658.8 hours (27.5 days)
- **Angular velocity**: 0.000152 deg/s

#### Comparative Analysis
- **Radius ratio (Moon/LEO)**: 58.4
- **Velocity ratio (LEO/Moon)**: 7.6
- **Period ratio (Moon/LEO)**: 446.7
- **LEO completes 1 orbit in 88.5 minutes**
- **Moon completes 0.0022 orbits in 88.5 minutes (0.81°)**

### Data Structure

The generated data arrays contain:
```python
leo_data = {
    'time_min': [0, 1, 2, 3, ...],           # Time in minutes
    'time_sec': [0, 60, 120, 180, ...],      # Time in seconds
    'theta_rad': [...],                       # Angular position in radians
    'theta_deg': [...],                       # Angular position in degrees
    'x_km': [...],                           # X coordinate in km
    'y_km': [...]                            # Y coordinate in km
}

moon_data = {
    # Same structure as leo_data
}
```

### Usage Examples

#### Access position at specific time:
```python
# Position at t = 60 minutes
time_index = 60
leo_angle = leo_data['theta_deg'][time_index]    # 244.1°
leo_x = leo_data['x_km'][time_index]             # -2,874 km
leo_y = leo_data['y_km'][time_index]             # -5,917 km
```

#### Calculate distance between objects:
```python
import math
distance = math.sqrt((moon_x - leo_x)**2 + (moon_y - leo_y)**2)
```

### Mathematical Foundation

The calculations are based on fundamental orbital mechanics:

1. **Circular orbital velocity**: v = √(μ/r)
2. **Orbital period**: T = 2π√(r³/μ)
3. **Angular velocity**: ω = √(μ/r³) = 2π/T
4. **Angular position**: θ(t) = ωt
5. **Cartesian coordinates**: 
   - x(t) = r·cos(θ(t))
   - y(t) = r·sin(θ(t))

Where:
- μ = 398,600.4418 km³/s² (Earth's gravitational parameter)
- r = orbital radius
- t = time from initial position

### File Locations
All files are located in:
```
lambert_solvers/duda_implementation/Porkchop/
├── porkchop.py           # Main calculator
├── detailed_analysis.py  # Mathematical analysis
├── data_access.py        # Data access examples
└── orbit_data.txt        # Complete dataset
```

### Features Implemented
✅ Circular orbit calculations for LEO satellite (200km altitude)  
✅ Circular orbit calculations for Moon (average Earth-Moon distance)  
✅ Coplanar, 2D orbital geometry  
✅ Position defined by angle θ with θ=0 at initial position  
✅ 1-minute time intervals  
✅ One complete LEO orbit duration (88.5 minutes)  
✅ Cartesian coordinate conversion  
✅ Data export functionality  
✅ Visualization capabilities  
✅ Mathematical verification  

The implementation provides a complete solution for defining and analyzing the position and time data for both circular orbits during one complete LEO satellite orbit as requested.