# Lunar Transfer Trajectory Analysis - Modular Structure

This project has been reorganized into a modular structure for better maintainability and code organization.

## Project Structure

```
Transfer_Simulation/
├── main.py                    # Main application entry point
├── config.py                  # Configuration and constants
├── interface.py               # User interface and I/O functions
├── earth_operations.py        # Earth departure calculations
├── trajectory_calculations.py # Geocentric trajectory analysis
├── lunar_operations.py        # Lunar SOI and orbit operations
├── analysis.py               # Parametric studies and optimization
├── plotting.py               # Visualization and plotting functions
└── Transfer.py               # Original monolithic file (legacy)
```

## Module Descriptions

### `main.py`
- **Purpose**: Main application entry point
- **Functions**: 
  - `main()`: Program entry point with menu system
  - `perform_single_calculation()`: Run single lambda1 analysis
  - `perform_parametric_study()`: Run parametric optimization study
- **Dependencies**: All other modules

### `config.py`
- **Purpose**: Central configuration and constants
- **Contents**: 
  - Physical constants (Earth/Moon radii, gravitational parameters)
  - Mission parameters (default values)
  - Conversion factors (DU to km, TU to seconds)
  - Analysis parameters (lambda1 range and step size)
- **Dependencies**: None

### `interface.py`
- **Purpose**: User interface and input/output functions
- **Functions**:
  - `get_user_input()`: Collect mission parameters from user
  - `display_organized_mission_summary()`: Show comprehensive mission results
  - `display_analysis_menu()`: Main menu system
- **Dependencies**: config.py

### `earth_operations.py`
- **Purpose**: Earth-related calculations
- **Functions**:
  - `calculate_earth_departure_delta_v()`: Departure maneuver analysis
- **Key Features**: 
  - Circular parking orbit analysis
  - Transfer trajectory properties
  - Delta-V calculations for Earth departure
- **Dependencies**: config.py

### `trajectory_calculations.py`
- **Purpose**: Geocentric trajectory analysis
- **Functions**:
  - `lunar_trajectory_calculations()`: Complete Earth-Moon transfer analysis
- **Key Features**:
  - Orbital mechanics calculations (energy, momentum, anomalies)
  - Time of flight analysis
  - Moon SOI intersection geometry
- **Dependencies**: config.py

### `lunar_operations.py`
- **Purpose**: Lunar sphere of influence and orbit operations
- **Functions**:
  - `lunar_soi_calculations()`: Lunar-centric trajectory analysis
  - `calculate_lunar_soi_transit_time()`: SOI entry to perigee time
  - `hyperbolic_to_elliptical_conversion()`: Orbit insertion maneuvers
- **Key Features**:
  - Hyperbolic trajectory analysis
  - Orbit conversion calculations
  - Circularization maneuvers
- **Dependencies**: config.py

### `analysis.py`
- **Purpose**: Parametric studies and optimization
- **Functions**:
  - `parametric_study_lambda1()`: Lambda1 optimization study
  - `find_optimal_lambda1()`: Find optimal solutions
  - `find_optimal_lambda1_complete()`: Complete mission optimization
- **Key Features**:
  - Batch calculations across lambda1 range
  - Multi-criteria optimization (delta-V, time, balanced)
  - Statistical analysis
- **Dependencies**: config.py, earth_operations.py, trajectory_calculations.py, lunar_operations.py

### `plotting.py`
- **Purpose**: Visualization and plotting
- **Functions**:
  - `create_parametric_plots()`: Standard parametric plots
  - `create_complete_parametric_plots()`: Comprehensive mission plots
- **Key Features**:
  - Multi-subplot visualizations
  - Optimization annotations
  - Publication-quality plots
- **Dependencies**: matplotlib, numpy

## Usage

### Running the Program
```bash
python main.py
```

### Menu Options
1. **Single Calculation**: Analyze one specific lambda1 value with detailed output
2. **Parametric Study**: Optimize lambda1 from 0-360° with plots and statistics

### Example Single Calculation
```python
from main import perform_single_calculation

# Run with user input
results = perform_single_calculation()
```

### Example Parametric Study
```python
from main import perform_parametric_study

# Run optimization study
lambda1_vals, delta_v_vals, time_vals, optimal = perform_parametric_study()
```

### Using Individual Modules
```python
from earth_operations import calculate_earth_departure_delta_v
from trajectory_calculations import lunar_trajectory_calculations
from lunar_operations import lunar_soi_calculations

# Calculate Earth departure
departure_results = calculate_earth_departure_delta_v(1.05, 1.372)

# Calculate geocentric trajectory
geo_results = lunar_trajectory_calculations(1.05, 1.372, 0, 30)

# Calculate lunar trajectory
lunar_results = lunar_soi_calculations(
    geo_results['r1'], geo_results['v1'], 
    geo_results['phi1_deg'], 30, geo_results['gamma1_deg']
)
```

## Configuration

Modify `config.py` to change:
- **Mission Parameters**: Parking orbit, transfer velocity, target altitude
- **Analysis Range**: Lambda1 range and step size
- **Physical Constants**: Gravitational parameters, celestial body properties

## Output Files

- `lambda1_parametric_study.png`: Standard parametric analysis plots
- `complete_mission_parametric_study.png`: Comprehensive mission analysis plots

## Migration from Legacy Code

The original `Transfer.py` file remains available for reference. The new modular structure provides:

- **Better Organization**: Related functions grouped by purpose
- **Easier Maintenance**: Changes isolated to relevant modules
- **Improved Reusability**: Individual modules can be imported separately
- **Cleaner Dependencies**: Clear separation of concerns
- **Enhanced Testability**: Individual components can be tested independently

## Dependencies

- `numpy`: Numerical calculations and array operations
- `matplotlib`: Plotting and visualization
- `math`: Mathematical functions

## Future Enhancements

The modular structure enables easy addition of:
- New optimization algorithms (add to `analysis.py`)
- Different plot types (add to `plotting.py`)
- Alternative mission profiles (modify `config.py`)
- Additional trajectory models (extend relevant modules)