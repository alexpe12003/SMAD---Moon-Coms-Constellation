# 🚀 Lambert Solvers

**A collection of Lambert's problem solvers and orbital transfer analysis tools for the SMAD Moon Communications Constellation project.**

Lambert's problem is a fundamental orbital mechanics problem that determines the orbit connecting two position vectors with a given time of flight. This is essential for calculating transfer trajectories between Earth and the Moon.

## 📁 Directory Structure

### 🎯 `duda_implementation/`
**Primary implementation directory containing modern Python-based Lambert solvers**

#### `Duda_Lambert/`
Core Lambert solver implementation with supporting utilities:
- **`lambert/lambert.py`** - Main Lambert solver using universal variable formulation
- **`Auxiliary/`** - Helper functions for coordinate transformations (COE ↔ State Vector)
- **`Test/`** - Test scripts and validation examples

#### `Porkchop/`
Comprehensive porkchop plot generator for Earth-Moon transfers:
- **`core_modules/`** - Essential orbit simulation and porkchop generation code
- **`final_results/`** - Generated porkchop plots and orbital data

### 📚 `Legacy/`
**Legacy implementations and alternative approaches**

#### `Main/`
- Earth-Moon transfer solver using Skyfield ephemeris data
- Solves Lambert's problem for real celestial body positions

## 🔧 Key Features

### Lambert Solver Capabilities
- **Universal Variable Formulation**: Robust numerical method that works for all orbit types
- **Prograde/Retrograde Options**: Support for both transfer directions
- **High Precision**: Configurable tolerance for convergence (default: 1e-8)
- **Multiple Orbit Types**: Handles elliptical, parabolic, and hyperbolic trajectories

### Porkchop Plot Analysis
- **LEO-to-Moon Transfers**: Analysis of Low Earth Orbit to lunar transfers
- **ΔV Optimization**: Total, departure, and arrival velocity requirements
- **Time Window Analysis**: Extended simulation periods (150+ hours)
- **Visual Output**: High-quality porkchop plots with optimal transfer identification

### Auxiliary Tools
- **Coordinate Conversions**: Classical Orbital Elements ↔ State Vectors
- **Orbit Propagation**: Position and velocity calculation over time
- **Real Ephemeris Data**: Integration with Skyfield for actual celestial positions

## 🚀 Quick Start

### Running the Main Lambert Solver
```python
# Navigate to the Duda_Lambert directory
cd duda_implementation/Duda_Lambert/Test/
python Test.py
```

### Generating Porkchop Plots
```python
# Navigate to the Porkchop core modules
cd duda_implementation/Porkchop/core_modules/
python quick_extended_porkchop.py
```

### Using Real Ephemeris Data
```python
# Navigate to the Legacy main directory
cd Legacy/Main/
python main.py
```

## 📊 Example Results

The porkchop analysis generates several types of plots:
- **Total ΔV Requirements**: Combined departure and arrival velocity needs
- **Departure ΔV**: Velocity change required to leave LEO
- **Arrival ΔV**: Velocity change required for lunar orbit insertion
- **Orbital Position Data**: Time-series data for LEO satellite and Moon positions

## 🔬 Technical Implementation

### Lambert Solver Algorithm
The core Lambert solver uses the universal variable formulation:
1. **Input**: Two position vectors (R₁, R₂) and time of flight (t)
2. **Process**: Newton's method iteration to solve for universal anomaly
3. **Output**: Initial and final velocity vectors (V₁, V₂)

### Key Equations
- Universal variable formulation for robust convergence
- Vis-viva equation for velocity calculations
- Stumpff functions for elliptical, parabolic, and hyperbolic cases

### Coordinate Systems
- **Position Vectors**: 3D Cartesian coordinates (km)
- **Velocity Vectors**: 3D Cartesian velocities (km/s)
- **Classical Elements**: Semi-major axis, eccentricity, inclination, etc.

## 🔧 Dependencies

### Python Packages
```
numpy          # Numerical computations
matplotlib     # Plotting and visualization
skyfield       # Astronomical calculations (Legacy implementation)
```

### Installation
```bash
pip install numpy matplotlib skyfield
```

## 🧪 Testing and Validation

The implementations include comprehensive test suites:
- **Analytical Validation**: Comparison with known orbital mechanics solutions
- **Convergence Testing**: Verification of numerical stability
- **Edge Case Handling**: Testing for extreme orbital parameters
- **Cross-Validation**: Comparison between different implementation approaches

## 📈 Applications

This Lambert solver collection is designed for:
- **Mission Planning**: Earth-Moon transfer trajectory design
- **Constellation Deployment**: Satellite positioning for communication networks
- **Orbital Mechanics Education**: Learning and research applications
- **Spacecraft Navigation**: Real-time trajectory calculations

## 🤝 Contributing

When working with these solvers:
1. Test thoroughly with known solutions
2. Validate convergence behavior
3. Document any modifications or improvements
4. Maintain compatibility with existing interfaces

## 📝 References

- Bate, R., Mueller, D., & White, J. (1971). *Fundamentals of Astrodynamics*
- Curtis, H. (2013). *Orbital Mechanics for Engineering Students*
- Vallado, D. (2013). *Fundamentals of Astrodynamics and Applications*

---

**Part of the SMAD Moon Communications Constellation Project**  
*Advanced orbital mechanics tools for lunar communication satellite deployment*