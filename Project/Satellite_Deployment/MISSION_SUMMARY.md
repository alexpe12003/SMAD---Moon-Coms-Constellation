# Moon Communications Constellation Deployment Summary

## Mission Overview
- **Starting orbit**: Elliptical 123km × 1500km at 62.8° inclination
- **Final constellation**: 42 satellites (7 carriers × 6 satellites each) in polar 1000km orbits
- **Total mission duration**: 28.6 hours (1.2 days)
- **Total mission Delta-V**: 8,745 m/s

## Mission Phases

### Phase 1: Orbit Circularization
- **Action**: Prograde burn at apogee to circularize orbit
- **Result**: 1500km circular orbit at 62.8° inclination
- **Delta-V**: 179.3 m/s (spacecraft only)
- **Time**: Performed at t = 0 (first apogee passage)

### Phase 2: Carrier Deployment
- **Action**: Deploy 7 carriers with equal time intervals
- **Orbit**: 1500km circular at 62.8° inclination
- **Interval**: 0.656 hours (2,361 seconds) between deployments
- **Duration**: 3.94 hours total
- **Delta-V**: 0 m/s (simple deployment)

### Phase 3: Inclination Changes
- **Action**: Each carrier changes inclination from 62.8° to 90° (polar)
- **Orbit**: Remains at 1500km circular altitude
- **Delta-V per carrier**: 578.7 m/s
- **Total Delta-V**: 4,051.2 m/s (7 carriers)
- **Time**: Performed immediately after each carrier deployment

### Phase 4: Satellite Deployment
- **Action**: Each carrier deploys 6 satellites to 1000km polar orbits
- **Transfer method**: Staggered Hohmann transfers (1500km → 1000km)
- **Staggering interval**: 2.674 hours between satellite deployments
- **Transfer time**: 2.035 hours per satellite
- **Delta-V per satellite**: 107.5 m/s
- **Total Delta-V**: 4,514.6 m/s (42 satellites)

## Delta-V Budget Summary

| Component | Unit Delta-V (m/s) | Quantity | Total Delta-V (m/s) |
|-----------|-------------------|----------|-------------------|
| **Spacecraft** | | | |
| Orbit circularization | 179.3 | 1 | 179.3 |
| **Carriers** | | | |
| Inclination change | 578.7 | 7 | 4,051.2 |
| **Satellites** | | | |
| Hohmann transfer | 107.5 | 42 | 4,514.6 |
| **TOTAL MISSION** | | | **8,745.0** |

## Timeline Summary

| Phase | Duration | Cumulative Time |
|-------|----------|----------------|
| Circularization | 0.0 hours | 0.0 hours |
| Carrier deployment | 3.94 hours | 3.94 hours |
| Inclination changes | Overlapped | 3.94 hours |
| Satellite deployment | 24.67 hours | 28.61 hours |
| **TOTAL MISSION** | **28.61 hours** | **(1.2 days)** |

## Final Constellation Configuration

- **Total satellites**: 42
- **Orbital planes**: 7 (polar orbits with different RAAN)
- **Satellites per plane**: 6
- **Satellite altitude**: 1000km
- **Satellite inclination**: 90° (polar)
- **Spacing within each plane**: 60°
- **Coverage**: Global lunar coverage
- **Revisit time**: ~0.59 hours per orbital plane

## Key Technical Details

### Initial Elliptical Orbit
- **Perigee**: 123km altitude (1,860.4km radius)
- **Apogee**: 1500km altitude (3,237.4km radius)
- **Semi-major axis**: 2,548.9km
- **Eccentricity**: 0.2701
- **Period**: 3.21 hours
- **Velocity at perigee**: 1.830 km/s
- **Velocity at apogee**: 1.051 km/s

### Circularized Orbit (1500km)
- **Velocity**: 1.231 km/s
- **Period**: 4.59 hours

### Final Satellite Orbit (1000km)
- **Velocity**: 1.633 km/s
- **Period**: 3.57 hours

### Satellite Transfer Parameters
- **Transfer semi-major axis**: 2,987.4km
- **First burn (retrograde at 1500km)**: 52.6 m/s
- **Second burn (prograde at 1000km)**: 54.9 m/s
- **Transfer time**: 2.035 hours

## Mission Advantages

1. **Efficient use of starting orbit**: Circularization at apogee maximizes efficiency
2. **Modular deployment**: Each carrier operates independently
3. **Optimal timing**: Staggered deployments achieve precise spacing
4. **Comprehensive coverage**: 7 polar planes provide global lunar coverage
5. **Reasonable Delta-V budget**: Total 8.7 km/s is achievable with modern propulsion
6. **Short mission duration**: Complete deployment in just over 1 day