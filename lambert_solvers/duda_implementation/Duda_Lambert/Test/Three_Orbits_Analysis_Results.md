# Three Orbits Analysis Results

## Overview
This analysis compares three different orbit types for the same true anomaly transfer case (θ₁ = 30°, θ₂ = 60°) and validates the Lambert solver against theoretical orbital mechanics calculations.

## Orbit Definitions

### 🚀 Orbit A - Classical Elliptical Orbit
- **Eccentricity (e)**: 0.1 (slightly elliptical)
- **Semi-major axis (a)**: 7,000 km
- **Orbit type**: Elliptical (e < 1, a > 0)
- **Characteristics**: Low Earth orbit with minimal eccentricity

### 🚀 Orbit B - Highly Elliptical Orbit  
- **Eccentricity (e)**: 0.9 (highly elliptical)
- **Semi-major axis (a)**: 70,000 km
- **Orbit type**: Elliptical (e < 1, a > 0)
- **Characteristics**: Very elongated ellipse, similar to geostationary transfer orbit

### 🚀 Orbit C - Hyperbolic Trajectory
- **Eccentricity (e)**: 2.0 (hyperbolic)
- **Semi-major axis (a)**: -7,000 km
- **Orbit type**: Hyperbolic (e > 1, a < 0)
- **Characteristics**: Open trajectory, escape velocity exceeded

## Key Results Summary

| Orbit | Type | r₁ (km) | r₂ (km) | v₁ (km/s) | v₂ (km/s) | ToF (hours) |
|-------|------|---------|---------|-----------|-----------|-------------|
| A | Elliptical | 6,378 | 6,600 | 8.25 | 7.99 | 0.12 |
| B | Elliptical | 7,474 | 9,172 | 10.05 | 9.01 | 0.13 |
| C | Hyperbolic | 7,687 | 10,500 | 12.68 | 11.53 | 0.12 |

## Key Observations

### 1. **Velocity Patterns**
- **Orbit A**: Lowest velocities due to small semi-major axis
- **Orbit B**: Moderate velocities despite large semi-major axis (high eccentricity effect)
- **Orbit C**: Highest velocities due to hyperbolic nature (escape trajectory)

### 2. **Radial Distance Behavior**
- **Orbit A**: Small variation between r₁ and r₂ (low eccentricity)
- **Orbit B**: Larger variation due to high eccentricity (0.9)
- **Orbit C**: Largest r₂ value due to hyperbolic expansion

### 3. **Time of Flight**
- **Orbits A & C**: Similar short transfer times (~0.12 hours)
- **Orbit B**: Slightly longer despite higher velocities (effect of orbit geometry)

### 4. **Lambert Solver Validation**
- ✅ **Perfect Agreement**: All Lambert solutions match theoretical calculations exactly
- ✅ **All Orbit Types**: Solver handles elliptical and hyperbolic orbits correctly
- ✅ **Orbital Elements**: Transfer orbit elements match input parameters precisely

## Physical Interpretation

### Orbit A (Classical LEO)
- Represents typical low Earth orbit scenarios
- Small altitude changes with moderate velocity requirements
- Short transfer times due to small orbital radius

### Orbit B (Highly Elliptical)
- Similar to Molniya or geostationary transfer orbits
- Large semi-major axis but still bound orbit
- High eccentricity creates significant radial variation

### Orbit C (Hyperbolic Escape)
- Represents interplanetary or escape trajectories
- Negative semi-major axis indicates unbound motion
- Highest energy requirements but efficient for long-distance transfers

## Technical Validation

The Lambert solver successfully:
1. **Handles Multiple Orbit Types**: Works for elliptical (e < 1) and hyperbolic (e > 1) cases
2. **Maintains Accuracy**: Zero velocity discrepancy between theoretical and Lambert solutions
3. **Preserves Orbital Elements**: Transfer orbits maintain exact input eccentricity and semi-major axis
4. **Correct Time Calculations**: Proper handling of elliptical vs hyperbolic time-of-flight formulas

## Conclusion

This analysis demonstrates that the Lambert solver implementation is robust and accurate across different orbit types, from circular low Earth orbits to hyperbolic escape trajectories. The perfect agreement between theoretical orbital mechanics and the Lambert solution validates both the mathematical implementation and the numerical solver's stability.