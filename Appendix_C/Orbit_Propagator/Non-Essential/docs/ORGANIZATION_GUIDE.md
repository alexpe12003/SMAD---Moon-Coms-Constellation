# Orbit Propagator - File Organization

## Core Implementation
- `Constants.py` - Physical constants and configurable parameters for Orbit B
- `Kepler.py` - Analytical orbit calculations and Kepler's equation solver
- `RK4.py` - Runge-Kutta 4th order numerical integration implementation  
- `RK8.py` - Runge-Kutta 8th order numerical integration implementation

## Main Programs
- `main_simulation.py` - Complete convergence study for multiple step sizes
- `rk_comparison.py` - Direct comparison between RK4 and RK8 methods

## Analysis Tools
- `analysis_utils.py` - Advanced analysis and plotting utilities

## Documentation
- `README.md` - Comprehensive user guide and documentation
- `project_summary.py` - Usage examples and project overview
- `ORGANIZATION_GUIDE.md` - This file (project structure)

## Results Directory
- `results/` - Auto-generated directory containing:
  - `convergence_summary.txt` - Tabulated convergence results
  - `orbit_analysis.png` - Comprehensive analysis plots
  - `rk4_results_step_Xs.npz` - Detailed numerical results for each step size
  - Additional plots from comparison studies

## Usage Patterns

### Quick Test
```bash
python Constants.py    # Test orbital parameters
python Kepler.py      # Test analytical solver
python RK4.py         # Test RK4 integrator
```

### Full Study
```bash
python main_simulation.py    # Complete convergence analysis
```

### Method Comparison  
```bash
python rk_comparison.py     # RK4 vs RK8 comparison
```

## Key Features
- Highly eccentric orbit (e=0.9) with challenging numerical properties
- Both 4th and 8th order Runge-Kutta implementations
- Comprehensive error analysis and convergence studies
- Energy and angular momentum conservation monitoring
- Flexible parameter configuration system
- Detailed plotting and visualization capabilities