# Lunar Transfer Trajectory Analysis - Organized Structure

This project has been reorganized into a modular folder structure for better maintainability and code organization.

## New Project Structure

```
Transfer_Simulation/
├── main.py                    # Main application entry point
├── core/                      # Core functionality
│   ├── __init__.py
│   ├── config.py             # Configuration and constants
│   └── interface.py          # User interface and I/O functions
├── operations/               # Calculation modules
│   ├── __init__.py
│   ├── earth_operations.py   # Earth departure calculations
│   ├── lunar_operations.py   # Lunar SOI and orbit operations
│   └── trajectory_calculations.py # Geocentric trajectory analysis
├── analysis/                 # Analysis and visualization
│   ├── __init__.py
│   ├── analysis.py          # Parametric studies and optimization
│   ├── plotting.py          # Visualization and plotting functions
│   ├── analyze_asin_error.py # Error analysis tools
│   └── diagnose_lambda1.py  # Lambda1 diagnostic tools
├── tests/                   # Test modules
│   ├── __init__.py
│   └── test_lunar_ops.py   # Unit tests for lunar operations
├── legacy/                  # Legacy code
│   └── Transfer.py         # Original monolithic file
└── docs/                   # Documentation
    ├── README.md           # This file
    └── Lambda1_Failure_Analysis.md
```

## Folder Descriptions

### `core/`
Contains the fundamental modules that provide configuration and interface functionality:
- **`config.py`**: Physical constants, mission parameters, conversion factors
- **`interface.py`**: User input/output functions and menu system

### `operations/`
Contains the core calculation modules for trajectory analysis:
- **`earth_operations.py`**: Earth departure delta-V calculations
- **`lunar_operations.py`**: Lunar SOI entry, transit time, and orbital mechanics
- **`trajectory_calculations.py`**: Geocentric trajectory analysis

### `analysis/`
Contains analysis, optimization, and visualization tools:
- **`analysis.py`**: Parametric studies and optimization algorithms
- **`plotting.py`**: Plotting and visualization functions
- **`analyze_asin_error.py`**: Diagnostic tools for mathematical domain errors
- **`diagnose_lambda1.py`**: Specific diagnostic tools for lambda1 parameter issues

### `tests/`
Contains unit tests and validation scripts:
- **`test_lunar_ops.py`**: Unit tests for lunar operations

### `legacy/`
Contains the original monolithic code for reference:
- **`Transfer.py`**: Original single-file implementation

### `docs/`
Contains project documentation and analysis reports:
- **`README.md`**: This documentation file
- **`Lambda1_Failure_Analysis.md`**: Detailed analysis of lambda1 parameter issues

## Usage

The main entry point remains the same:

```bash
python main.py
```

The modular structure allows for easier:
- **Development**: Each module has a specific responsibility
- **Testing**: Individual modules can be tested in isolation
- **Maintenance**: Bug fixes and improvements are localized
- **Extension**: New features can be added without affecting existing code

## Import Structure

The new modular structure uses relative imports to maintain clean dependencies:
- Core modules are imported from `core/`
- Operations are imported from `operations/`
- Analysis tools are imported from `analysis/`

This organization makes the codebase more professional and easier to navigate while maintaining all existing functionality.