# Lunar Transfer Trajectory Analysis

This module provides a comprehensive lunar transfer trajectory analysis system for calculating Earth-Moon transfer trajectories with various optimization capabilities.

## Overview

The Transfer_Simulation module is designed to analyze lunar transfer trajectories using various parameters and optimization techniques. It provides both single calculations and parametric studies to find optimal mission parameters for lunar transfers.

## Features

- **Single Trajectory Calculation**: Calculate a single lunar transfer with user-defined parameters
- **Parametric Studies**: Analyze lambda1 parameter variations from 0-360 degrees
- **Optimal Parameters Calculation**: Run calculations with pre-optimized parameters
- **Multi-Parameter Optimization**: Optimize multiple parameters to minimize mission delta-V
- **Fixed V0 Optimization**: Optimization with fixed velocity parameters

## Structure

```
Transfer_Simulation/
├── main.py                     # Main application entry point
├── core/                       # Core interface modules
│   ├── __init__.py
│   ├── interface.py            # User interface and menu system
│   └── config.py               # Configuration parameters
├── operations/                 # Mission operations calculations
│   ├── __init__.py
│   ├── earth_operations.py     # Earth departure calculations
│   ├── trajectory_calculations.py  # Trajectory computations
│   └── lunar_operations.py     # Lunar arrival and orbit calculations
└── analysis/                   # Analysis and optimization tools
    ├── __init__.py
    ├── analysis.py             # Parametric studies and analysis
    ├── plotting.py             # Visualization and plotting
    └── optimization.py         # Multi-parameter optimization
```

## Usage

### Running the Main Application

```bash
python main.py
```

The application will present a menu with the following options:

1. **Single Calculation**: Perform a single lunar transfer calculation with custom parameters
2. **Parametric Study**: Run a parametric study of lambda1 values (0-360°)
3. **Optimal Parameters**: Use pre-calculated optimal parameters for the mission
4. **Multi-Parameter Optimization**: Optimize V0, gamma0, and lambda1 parameters
5. **Fixed V0 Optimization**: Optimize with V0 fixed at 1.372 DU/TU

### Input Parameters

- **R0**: Initial orbital radius (DU - Distance Units)
- **V0**: Initial velocity (DU/TU - Distance Units per Time Unit)
- **gamma0**: Flight path angle at departure (degrees)
- **lambda1**: True anomaly at lunar SOI entry (degrees)

### Output Results

The system provides comprehensive mission analysis including:

- Earth departure delta-V requirements
- Geocentric trajectory parameters
- Lunar SOI entry conditions
- Lunar orbit insertion parameters
- Mission duration and timing
- Visualization plots (for parametric studies)

## Key Functions

### Main Application Functions

- `perform_single_calculation()`: Single trajectory analysis
- `perform_parametric_study()`: Lambda1 parametric study with plotting
- `perform_optimal_parameters_calculation()`: Uses optimized parameters
- `perform_multi_parameter_optimization()`: Full parameter optimization
- `perform_fixed_v0_optimization()`: Constrained optimization

### Core Modules

- **Interface**: User input/output, menu system, result display
- **Earth Operations**: Launch and departure calculations
- **Trajectory Calculations**: Geocentric transfer computations
- **Lunar Operations**: Lunar SOI and orbit insertion calculations
- **Analysis**: Parametric studies and optimal parameter finding
- **Plotting**: Visualization of results and parametric studies
- **Optimization**: Multi-parameter optimization algorithms

## Dependencies

The system requires the following Python packages:
- `numpy` - Numerical computations
- `matplotlib` - Plotting and visualization
- `scipy` - Optimization algorithms
- Standard Python libraries: `sys`, `os`, `math`

## Coordinate System and Units

- **Distance Units (DU)**: Earth radii (1 DU = 6,378 km)
- **Time Units (TU)**: Canonical time units
- **Velocity Units**: DU/TU
- **Angles**: Degrees for user interface, radians for calculations

## Optimal Parameters

The system includes pre-calculated optimal parameters from optimization runs:

- **R0**: 1.050 DU
- **V0**: 1.372 DU/TU
- **γ0**: 0.0°
- **λ1**: 236.4°

These parameters provide collision-free trajectories with minimized delta-V requirements.

## Example Workflow

1. Run `python main.py`
2. Select option from the menu
3. For single calculations: Enter R0, V0, gamma0, and lambda1 values
4. For optimization: Specify number of optimization steps
5. Review comprehensive mission summary with all trajectory parameters
6. For parametric studies: View generated plots showing parameter relationships

## Notes

- All calculations use canonical units for computational efficiency
- Results include both dimensional and non-dimensional values
- The system performs collision detection and warns of potential issues
- Optimization algorithms use gradient-free methods suitable for aerospace applications
- Visualization capabilities help understand parameter sensitivities and trade-offs

## Contributing

When modifying the code:
1. Maintain the modular structure
2. Update unit tests if available
3. Ensure compatibility with the main application interface
4. Document any new parameters or calculation methods