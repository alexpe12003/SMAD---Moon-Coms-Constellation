# Orbit B Numerical Integration Study

This project implements a comprehensive numerical integration study for satellite orbital mechanics using Runge-Kutta methods. The study focuses on **Orbit B** with an eccentricity of 0.9 and a semi-major axis of 70,000 km.

## Project Structure

```
Orbit_Propagator/
├── Constants.py           # Physical constants and configurable parameters
├── Kepler.py             # Analytical orbit calculations and Kepler's equation solver
├── RK4.py                # Runge-Kutta 4th order numerical integration
├── RK8.py                # Runge-Kutta 8th order numerical integration (optional)
├── main_simulation.py    # Main simulation script
├── analysis_utils.py     # Analysis and plotting utilities
├── README.md            # This file
└── results/             # Output directory (created automatically)
```

## Features

### 1. Orbit B Parameters
- **Semi-major axis**: 70,000 km
- **Eccentricity**: 0.9 (highly eccentric orbit)
- **Inclination**: 0° (equatorial orbit)
- **Periapsis distance**: 7,000 km
- **Apoapsis distance**: 133,000 km
- **Orbital period**: ~22.7 hours

### 2. Numerical Integration
- **RK4 Implementation**: 4th order Runge-Kutta with fixed and adaptive step sizes
- **RK8 Implementation**: 8th order Runge-Kutta for high-precision applications
- **Multiple Step Sizes**: Tests convergence with steps from 1s to 300s
- **Conservation Analysis**: Monitors energy and angular momentum conservation
- **Error Analysis**: Compares numerical solution with analytical reference

### 3. Robust Orbital Elements Calculation
- **Integrated coe_sv.py**: Uses Algorithm 4.1 from orbital mechanics literature
- **Special Cases Handling**: Proper treatment of circular and equatorial orbits
- **High Numerical Stability**: Robust quadrant determination and error checking
- **Comprehensive Testing**: Validated across all orbital regimes

### 3. Analysis Capabilities
- Convergence order estimation
- Error metrics calculation
- Conservation property analysis
- Computational efficiency analysis
- Comprehensive plotting and visualization
- Results export to CSV and markdown reports

## Quick Start

### Running the Simulation

1. **Basic simulation run:**
```bash
cd Orbit_Propagator
python main_simulation.py
```

2. **Test individual components:**
```bash
# Test constants and orbital parameters
python Constants.py

# Test Kepler equation solver
python Kepler.py

# Test RK4 integrator
python RK4.py

# Test analysis utilities
python analysis_utils.py
```

### Configuration

All configurable parameters are in `Constants.py`:

- **Orbit parameters**: Modify orbital elements in the constants section
- **Integration parameters**: Adjust step sizes and study duration
- **Analysis parameters**: Configure tolerances and output options

```python
# Example modifications in Constants.py
ORBIT_B_ECCENTRICITY = 0.5        # Change eccentricity
INTEGRATION_STEPS = [10, 60, 300]  # Use different step sizes
NUM_ORBITS_STUDY = 3              # Integrate for 3 orbits instead of 5
```

## Results and Output

The simulation generates:

1. **Console Output**: Real-time progress and summary statistics
2. **Results Directory**: 
   - `convergence_summary.txt`: Tabulated convergence results
   - `rk4_results_step_Xs.npz`: Detailed numerical results for each step size
   - `orbit_analysis.png`: Comprehensive analysis plots
   - `convergence_results.csv`: Exportable data
   - `analysis_report.md`: Detailed markdown report

3. **Plots**:
   - 3D orbital trajectories
   - Convergence analysis (log-log plots)
   - Energy conservation monitoring
   - Computational efficiency analysis

## Understanding the Results

### Convergence Analysis
The simulation tests different step sizes and analyzes how errors decrease as step size decreases. For RK4, you should expect:
- **Convergence order**: ~4.0 (theoretical)
- **Error behavior**: Error ∝ (step size)^4

### Energy Conservation
The two-body problem conserves energy. Monitor energy errors to assess:
- **Numerical stability**: Energy errors should remain bounded
- **Long-term behavior**: Energy drift indicates integration issues
- **Step size effects**: Smaller steps generally improve conservation

### Optimal Step Size
The analysis identifies optimal step sizes based on:
- **Accuracy requirements**: Maximum acceptable error
- **Computational budget**: Available computation time
- **Efficiency metric**: Error per unit computation time

## Advanced Usage

### Custom Orbital Elements
To study different orbits, modify `Constants.py`:

```python
# Example: Circular orbit at 400 km altitude
ORBIT_B_SEMIMAJOR_AXIS = 6371 + 400  # km
ORBIT_B_ECCENTRICITY = 0.0            # circular

# Example: Moon transfer orbit
ORBIT_B_SEMIMAJOR_AXIS = 200000.0     # km
ORBIT_B_ECCENTRICITY = 0.97           # highly eccentric
```

### Adding Perturbations
The RK4 implementation supports perturbations. Example usage:

```python
def atmospheric_drag(t, state):
    """Simple atmospheric drag model."""
    r = state[:3]
    v = state[3:]
    # Implement drag acceleration
    return drag_acceleration

# Use in propagation
times, pos, vel = propagate_orbit_rk4(
    initial_state, time_span, step_size,
    perturbations=[atmospheric_drag]
)
```

### Custom Analysis
Use `analysis_utils.py` functions for custom analysis:

```python
from analysis_utils import calculate_position_error_metrics, estimate_convergence_order

# Analyze custom results
metrics = calculate_position_error_metrics(numerical_pos, analytical_pos, times)
order, r_squared, _ = estimate_convergence_order(step_sizes, errors)
```

## Theoretical Background

### Runge-Kutta Methods
The RK4 method approximates the solution to:
```
dy/dt = f(t, y)
```

Using the update formula:
```
y_{n+1} = y_n + (k1 + 2k2 + 2k3 + k4)/6
```

Where:
- k1 = h*f(t_n, y_n)
- k2 = h*f(t_n + h/2, y_n + k1/2)
- k3 = h*f(t_n + h/2, y_n + k2/2)
- k4 = h*f(t_n + h, y_n + k3)

### Two-Body Problem
The equations of motion for the two-body problem:
```
d²r/dt² = -μr/|r|³
```

Where:
- r: position vector
- μ: gravitational parameter (GM)

### Kepler's Equation
The analytical solution uses Kepler's equation:
```
M = E - e*sin(E)
```

Where:
- M: mean anomaly
- E: eccentric anomaly  
- e: eccentricity

## Validation and Verification

The implementation includes several validation checks:

1. **Kepler Equation Solver**: Verified against known solutions
2. **RK4 Integrator**: Tested with harmonic oscillator (analytical solution)
3. **Conservation Laws**: Energy and angular momentum monitoring
4. **Analytical Comparison**: Direct comparison with analytical orbit propagation

## Performance Considerations

- **Step Size Selection**: Balance between accuracy and computation time
- **Adaptive Integration**: RK4 adaptive step size control available
- **Memory Usage**: Large time spans require careful memory management
- **Numerical Precision**: Double precision recommended for high eccentricity orbits

## Troubleshooting

### Common Issues

1. **Import Errors**: Ensure all files are in the same directory
2. **Memory Issues**: Reduce `NUM_ORBITS_STUDY` or increase step sizes
3. **Convergence Problems**: Check orbital parameters for physical validity
4. **Plot Display**: Matplotlib backend issues may prevent plot display

### Error Messages

- **"Kepler equation did not converge"**: Increase tolerance or max iterations
- **"Invalid orbital elements"**: Check for negative semi-major axis or e > 1
- **"Memory error"**: Reduce integration time span or increase step size

## Extensions and Modifications

### Potential Enhancements
1. **Higher-Order Methods**: Implement RK8 or adaptive methods
2. **Perturbation Models**: Add atmospheric drag, J2, solar radiation pressure
3. **Multi-Body Problem**: Extend to restricted three-body problem
4. **Parallel Processing**: Parallelize convergence study across step sizes

### Research Applications
- Satellite mission design
- Orbital mechanics education
- Numerical methods comparison
- Spacecraft trajectory optimization

## References

1. Curtis, H.D. "Orbital Mechanics for Engineering Students"
2. Vallado, D.A. "Fundamentals of Astrodynamics and Applications"
3. Battin, R.H. "An Introduction to the Mathematics and Methods of Astrodynamics"

## License

This code is generated for educational purposes as part of the SMAD Moon Communications Constellation Study.