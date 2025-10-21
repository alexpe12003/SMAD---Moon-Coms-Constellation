# New Features Added to Orbit Simulation

## Overview
I have successfully added the requested functionality to the Orbit B numerical integration study. The enhanced simulation now provides comprehensive options for running individual methods (RK4 or RK8) with detailed plotting capabilities.

## New Menu Options

### 1. Enhanced Interactive Menu
The main menu now offers 6 options instead of 4:

1. **Run Full Convergence Study (Both Methods)** - Original functionality
2. **Run RK4 Method Only** - NEW
3. **Run RK8 Method Only** - NEW  
4. **Run Special Case Only** - Original functionality
5. **Help & Information** - Original functionality
6. **Exit** - Original functionality

### 2. New Command Line Options
In addition to the interactive menu, you can now run specific methods directly:

```bash
python main_simulation.py --rk4    # Run RK4 method only
python main_simulation.py --rk8    # Run RK8 method only
python main_simulation.py --full   # Run both methods (original)
python main_simulation.py --special-case  # Run special case only
```

## New Plotting Features

For each individual method (RK4 or RK8), the simulation now creates comprehensive plots:

### Position, Velocity, and True Anomaly vs Time Plots
- **Position Magnitude vs Time**: Shows analytical vs numerical solutions for 4 different step sizes (1s, 30s, 100s, 300s)
- **Velocity Magnitude vs Time**: Shows analytical vs numerical solutions for 4 different step sizes  
- **True Anomaly vs Time**: Shows analytical vs numerical solutions for 4 different step sizes

Each plot includes:
- Black dashed line for analytical (reference) solution
- Colored lines for numerical solutions at different step sizes
- Error statistics displayed as text boxes (max error, final error)
- Time axis in hours for readability

### Error Analysis Plots
- **Position Errors vs Time**: Semi-log plot showing position errors over time for ALL step sizes
- **Velocity Errors vs Time**: Semi-log plot showing velocity errors over time for ALL step sizes  
- **True Anomaly Errors vs Time**: Semi-log plot showing true anomaly errors over time for ALL step sizes

Each error plot shows:
- Different colored lines for each step size
- Maximum error values in the legend
- Time evolution of numerical errors
- Semi-logarithmic scale for clear error visualization

### Additional Analysis Plots
- **Convergence Analysis**: Log-log plots showing error vs step size with theoretical convergence lines
- **3D Trajectory Comparison**: 3D visualization comparing analytical and numerical trajectories
- **Computation Time vs Accuracy Trade-off**: Shows performance vs accuracy relationship

## File Structure

### Generated Output Files
When running individual methods, the following files are created in the `results/` directory:

**For RK4 Method:**
- `rk4_method_summary.txt` - Text summary of results
- `rk4_single_method_step_[X]s.npz` - Numerical data for each step size
- `rk4_position_vs_time.png` - Position vs time plots
- `rk4_velocity_vs_time.png` - Velocity vs time plots  
- `rk4_true_anomaly_vs_time.png` - True anomaly vs time plots
- `rk4_position_errors_vs_time.png` - Position error analysis
- `rk4_velocity_errors_vs_time.png` - Velocity error analysis
- `rk4_true_anomaly_errors_vs_time.png` - True anomaly error analysis
- `rk4_convergence_analysis.png` - Convergence study plots
- `rk4_3d_trajectory.png` - 3D trajectory comparison

**For RK8 Method:**
- Same file structure as RK4 but with `rk8_` prefix

### Data Storage Format
Each `.npz` file contains:
- `times` - Time points array
- `positions` - Numerical position vectors [N×3]
- `velocities` - Numerical velocity vectors [N×3] 
- `analytical_positions` - Analytical position vectors [N×3]
- `analytical_velocities` - Analytical velocity vectors [N×3]
- `analytical_true_anomaly` - Analytical true anomaly values
- `position_errors` - Position error magnitudes
- `velocity_errors` - Velocity error magnitudes  
- `true_anomaly_method` - Numerical true anomaly values
- `true_anomaly_errors` - True anomaly error magnitudes
- `step_size` - Integration step size used
- `method` - Integration method name

## Technical Implementation Details

### New Class Methods Added
To the `OrbitSimulation` class:
- `run_single_method_convergence_study(method_name)` - Runs convergence study for one method
- `analyze_single_method_results(method_name)` - Analyzes and prints results for one method  
- `save_single_method_results(method_name)` - Saves results to files for one method

### New Plotting Functions Added
To `plotting_utils.py`:
- `create_single_method_plots()` - Main plotting coordinator for individual methods
- `create_position_vs_time_plots()` - Position magnitude vs time comparison plots
- `create_velocity_vs_time_plots()` - Velocity magnitude vs time comparison plots
- `create_true_anomaly_vs_time_plots()` - True anomaly vs time comparison plots
- `create_position_error_vs_time_plots()` - Position error evolution plots
- `create_velocity_error_vs_time_plots()` - Velocity error evolution plots
- `create_true_anomaly_error_vs_time_plots()` - True anomaly error evolution plots
- `create_single_method_convergence_plot()` - Convergence analysis plots
- `create_single_method_3d_trajectory()` - 3D trajectory comparison plots

## Usage Examples

### Running RK4 Method Only
```bash
# Command line
python main_simulation.py --rk4

# Or interactive menu: select option 2
python main_simulation.py
```

This will:
1. Test RK4 with step sizes: 1, 5, 10, 30, 60, 100, 300 seconds
2. Run simulation for 5 complete orbits (~10.6 days)
3. Generate 8 comprehensive plot files showing analytical vs numerical comparisons
4. Save detailed numerical results for each step size
5. Provide convergence analysis and error statistics

### Running RK8 Method Only  
```bash
# Command line
python main_simulation.py --rk8

# Or interactive menu: select option 3
python main_simulation.py
```

Same functionality as RK4 but using the 8th-order Runge-Kutta method.

## Key Benefits

1. **Focused Analysis**: You can now analyze just one method without running both
2. **Comprehensive Visualization**: Each method gets 8 detailed plots showing all aspects of the solution
3. **Error Evolution**: See how errors grow over time for each step size
4. **Quick Comparison**: Side-by-side analytical vs numerical solutions at multiple step sizes
5. **Performance Analysis**: Understand computation time vs accuracy trade-offs
6. **Detailed Data**: All numerical results saved for further analysis

## Backward Compatibility

All original functionality remains unchanged:
- Option 1 (Full Convergence Study) works exactly as before
- Option 4 (Special Case) works exactly as before  
- All original command line options still work
- Original file outputs are preserved

The new features are purely additive and don't affect existing workflows.