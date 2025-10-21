# Multi-Parameter Optimization for Lunar Transfer Missions

## Overview

This module provides a comprehensive optimization function that varies multiple mission parameters to find the optimal combination that minimizes the total mission delta-V for lunar transfer trajectories.

## Features

- **Multi-parameter optimization**: Simultaneously optimizes V0 (initial velocity), gamma0 (flight path angle), and lambda1 (lunar SOI crossing angle)
- **Fixed R0**: Keeps the parking orbit radius fixed at 1.05 DU as specified
- **Comprehensive analysis**: Includes complete mission delta-V calculation (Earth departure + lunar operations)
- **Progress tracking**: Real-time progress updates and success rate monitoring
- **Display only**: Results are displayed on screen without creating files
- **Statistical analysis**: Provides comprehensive statistics on the optimization results

## Parameter Ranges

| Parameter | Range | Units | Description |
|-----------|-------|--------|-------------|
| R0 | 1.05 (fixed) | DU | Parking orbit radius |
| V0 | 1.372 - 2.0 | DU/TU | Initial transfer velocity |
| gamma0 | 0 - 90 | degrees | Flight path angle at departure |
| lambda1 | 0 - 360 | degrees | Lunar SOI crossing angle |

## Usage

### Method 1: Through Main Interface

```python
# Run the main program and select option 3
python main.py
```

### Method 2: Direct Optimization

```python
# Run the demonstration with exactly 360 steps
python demo_optimization.py
```

### Method 3: Custom Optimization

```python
from analysis.optimization import multi_parameter_optimization

# Run optimization with custom parameters
results = multi_parameter_optimization(num_steps=500)

# Results are displayed automatically (no file saving)
```

## Function Details

### `multi_parameter_optimization(num_steps=360)`

**Parameters:**
- `num_steps` (int): Minimum number of optimization steps. The function uses the cube root to distribute steps across the 3D parameter space.

**Returns:**
- Dictionary containing:
  - `optimal_result`: Best parameter combination found
  - `all_results`: Complete list of all successful calculations
  - `statistics`: Optimization statistics and success rates

**Example Output:**
```
OPTIMAL MISSION PARAMETERS
======================================================================
Minimum Mission Delta-V: 4.835 km/s
Optimal Parameters:
  R0 = 1.050 DU
  V0 = 1.372 DU/TU
  γ0 = 0.0°
  λ1 = 0.0°

Mission Breakdown:
  Earth Departure ΔV: 3.131 km/s
  Lunar Operations ΔV: 1.703 km/s
  Total Mission ΔV: 4.835 km/s
  Total Mission Time: 59.7 hours
```

**Note**: File saving functionality has been removed. Results are displayed on screen only.

## Performance Notes

- **Calculation Time**: Each parameter combination requires multiple trajectory calculations. Expect 1-5 seconds per successful calculation.
- **Success Rate**: Typically 50-80% of parameter combinations result in valid trajectories.
- **Memory Usage**: Results are stored in memory. For large optimizations (>10,000 steps), consider periodic saving.

## Optimization Strategy

The function uses a grid search approach:
1. Distributes the specified number of steps across the 3D parameter space
2. Creates uniform grids for V0, gamma0, and lambda1
3. Evaluates every combination in the grid
4. Tracks the best result in real-time
5. Provides comprehensive analysis of the parameter space

## Example Results

For a 360-step optimization, you can expect to find:
- **Optimal V0**: Typically around 1.372-1.4 DU/TU (minimum energy transfers)
- **Optimal gamma0**: Often 0° (horizontal departure)
- **Optimal lambda1**: Varies significantly, typically 0-30° or 330-360°
- **Delta-V savings**: 0.1-0.5 km/s compared to non-optimized parameters

## File Generation

No files are created during optimization. All results are displayed on screen only.

## Integration with Main Program

The optimization function is fully integrated into the main lunar transfer analysis program:

1. Run `main.py`
2. Select option "3. Multi-parameter optimization: V0, gamma0, lambda1"
3. Enter the number of optimization steps (minimum 360)
4. Wait for optimization to complete
5. View optimal parameters and mission breakdown

This provides a comprehensive tool for finding the most efficient lunar transfer trajectory parameters within the specified constraints.