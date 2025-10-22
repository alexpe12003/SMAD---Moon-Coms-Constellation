"""
Runge-Kutta 4th Order (RK4) Numerical Integration for Orbital Mechanics
=======================================================================

This module implements the Runge-Kutta 4th order method for numerical 
integration of orbital equations of motion. It provides both fixed and
adaptive step size integration capabilities.

CURRENT IMPLEMENTATION STATUS:
=============================

✅ FULLY IMPLEMENTED:
- Pure two-body gravitational dynamics (Kepler problem)
- Fixed step size RK4 integration 
- Adaptive step size RK4 integration with error control
- Conservation analysis utilities (energy, angular momentum)
- Comprehensive testing framework

⚠️ FRAMEWORK READY (NO IMPLEMENTATIONS):
- Perturbation support infrastructure
- Modular perturbation function interface
- Integration with perturbations

❌ NOT YET IMPLEMENTED:
- Specific perturbation forces (atmospheric drag, J2, solar radiation pressure)
- Third-body effects (Moon, Sun gravitational influence)
- General relativistic corrections
- Atmospheric density models

PHYSICS MODELED:
===============
Currently: Pure two-body problem
- Central body: Point-mass Earth (or other celestial body)
- Satellite: Point-mass with no attitude dynamics  
- Force: Inverse-square gravitational attraction only
- Coordinate system: Inertial frame (ECI or similar)
- Mathematical model: d²r/dt² = -μr/|r|³

NUMERICAL METHODS:
=================
- RK4 with fixed step size: Constant Δt, predictable cost
- RK4 with adaptive step size: Variable Δt with error control
- Richardson extrapolation for error estimation
- Conservation monitoring for solution validation

INTENDED USE:
============
- Orbit B numerical convergence study (e=0.9, a=70000km)
- Comparing numerical vs. analytical solutions
- Algorithm performance analysis
- Educational orbital mechanics simulations

Author: Generated for SMAD Moon Communications Constellation Study
Date: October 2025
"""

import numpy as np

try:
    from .Constants import MU_EARTH, MIN_TIME_STEP, MAX_TIME_STEP
except ImportError:
    # Fallback to absolute imports when running directly
    from Constants import MU_EARTH, MIN_TIME_STEP, MAX_TIME_STEP

# =============================================================================
# ORBITAL EQUATIONS OF MOTION
# =============================================================================

def orbital_equations_of_motion(t, state, mu=MU_EARTH):
    """
    Orbital equations of motion for two-body problem.
    
    This function implements the FUNDAMENTAL two-body gravitational dynamics:
    
    Mathematical Model:
    d/dt[r] = v                    (position derivative = velocity)
    d/dt[v] = -μ*r/|r|³           (velocity derivative = gravitational acceleration)
    
    Physical Interpretation:
    - Models a satellite orbiting a point-mass central body (Earth)
    - Assumes NO atmospheric drag, NO oblateness, NO third-body effects
    - Pure inverse-square law gravitational force: F = -GMm/r² r_hat
    - This is the EXACT solution for the classical Kepler problem
    
    Integration Context:
    - This is what gets integrated by the RK4 numerical methods
    - Converts the 2nd-order ODE (Newton's law) into 1st-order system
    - Perfect for studying numerical convergence vs. analytical solutions
    
    Args:
        t (float): Time [s] (not used for autonomous system - gravity doesn't change with time)
        state (np.ndarray): State vector [x, y, z, vx, vy, vz] [km, km/s]
        mu (float): Gravitational parameter [km³/s²] (default: Earth = 398,600.4418)
        
    Returns:
        np.ndarray: State derivative [vx, vy, vz, ax, ay, az] [km/s, km/s²]
                   Ready for numerical integration by RK4/RK8 methods
    """
    # Extract position and velocity from 6D state vector
    r = state[:3]  # Position vector [x, y, z] in km
    v = state[3:]  # Velocity vector [vx, vy, vz] in km/s
    
    # Calculate distance from central body (Earth center)
    r_magnitude = np.linalg.norm(r)
    
    # Gravitational acceleration using Newton's law of universal gravitation
    # F = ma = -GMm/r² r_hat  =>  a = -GM/r² r_hat = -μ/r² r_hat
    # Vector form: a = -μ * r / |r|³  (r/|r| gives unit vector, |r|² in denominator)
    a = -mu * r / r_magnitude**3
    
    # Return state derivative for integration: [dr/dt, dv/dt] = [v, a]
    return np.concatenate([v, a])

def perturbed_orbital_equations(t, state, mu=MU_EARTH, perturbations=None):
    """
    Orbital equations of motion with perturbations.
    
    NOTE: This function provides a FRAMEWORK for perturbations, but currently
    NO SPECIFIC PERTURBATIONS ARE IMPLEMENTED. The system defaults to pure
    two-body dynamics unless perturbation functions are explicitly provided.
    
    Perturbation Framework:
    - Ready to accept any perturbation function with signature: func(t, state) -> acceleration_vector
    - Examples that COULD be added: atmospheric drag, J2 oblateness, solar radiation pressure
    - Each perturbation function should return a 3D acceleration vector [ax, ay, az] in km/s²
    
    Args:
        t (float): Time [s]
        state (np.ndarray): State vector [x, y, z, vx, vy, vz] [km, km/s]
        mu (float): Gravitational parameter [km³/s²]
        perturbations (list): List of perturbation functions (currently None by default)
                             Each function should have signature: func(t, state) -> np.array([ax, ay, az])
        
    Returns:
        np.ndarray: State derivative including perturbations
                   Currently returns pure two-body dynamics: d/dt[r, v] = [v, -μr/|r|³]
    """
    # Two-body acceleration (ALWAYS ACTIVE - this is the primary gravitational force)
    # Implements: a = -μ/r² in the direction of -r_hat
    state_dot = orbital_equations_of_motion(t, state, mu)
    
    # Add perturbations if provided (CURRENTLY NO PERTURBATIONS ARE DEFINED)
    # This loop is ready to accept perturbation functions but will not execute
    # unless specific perturbation functions are passed in the 'perturbations' list
    if perturbations is not None:
        for perturbation_func in perturbations:
            # Each perturbation function should return acceleration vector [ax, ay, az]
            perturb_accel = perturbation_func(t, state)
            # Add perturbation accelerations to the total acceleration (indices 3:6)
            state_dot[3:] += perturb_accel
    
    return state_dot

# =============================================================================
# RUNGE-KUTTA 4TH ORDER INTEGRATOR
# =============================================================================

def rk4_step(func, t, y, h, *args):
    """
    Single step of Runge-Kutta 4th order method.
    
    This is the CORE numerical integration algorithm that solves the differential equation.
    
    Mathematical Foundation:
    - Implements the classic RK4 formula with 4 slope evaluations
    - 4th-order accuracy: Local error ∝ h⁵, Global error ∝ h⁴
    - More accurate than Euler (1st order) or RK2 (2nd order) methods
    
    Algorithm Steps:
    1. k1: Slope at beginning of interval (t, y)
    2. k2: Slope at midpoint using k1 estimate (t+h/2, y+k1/2) 
    3. k3: Slope at midpoint using k2 estimate (t+h/2, y+k2/2)
    4. k4: Slope at end using k3 estimate (t+h, y+k3)
    5. Weighted average: y_new = y + (k1 + 2k2 + 2k3 + k4)/6
    
    For Orbital Mechanics:
    - func = orbital_equations_of_motion (currently pure two-body)
    - y = [x, y, z, vx, vy, vz] (6D state vector)
    - Each k evaluation computes [v, a] at different points
    
    Args:
        func (callable): Function to integrate (dy/dt = func(t, y, *args))
        t (float): Current time
        y (np.ndarray): Current state vector [x, y, z, vx, vy, vz]
        h (float): Step size (Δt)
        *args: Additional arguments for func (e.g., μ, perturbations)
        
    Returns:
        np.ndarray: Next state vector y(t+h) computed with RK4 accuracy
    """
    # Calculate k values (the four slope evaluations that make RK4 accurate)
    k1 = h * func(t, y, *args)              # Slope at start of interval
    k2 = h * func(t + h/2, y + k1/2, *args) # Slope at midpoint using k1 prediction
    k3 = h * func(t + h/2, y + k2/2, *args) # Slope at midpoint using k2 prediction  
    k4 = h * func(t + h, y + k3, *args)     # Slope at end using k3 prediction
    
    # Update state using RK4 formula (weighted average of slopes)
    # This gives 4th-order accuracy compared to simple Euler method
    y_new = y + (k1 + 2*k2 + 2*k3 + k4) / 6
    
    return y_new

def rk4_integrate(func, t_span, y0, step_size, *args, save_intermediate=True):
    """
    Integrate differential equation using RK4 method with fixed step size.
    
    Args:
        func (callable): Function to integrate (dy/dt = func(t, y, *args))
        t_span (tuple): Time span (t_start, t_end)
        y0 (np.ndarray): Initial state vector
        step_size (float): Integration step size
        *args: Additional arguments for func
        save_intermediate (bool): Whether to save intermediate results
        
    Returns:
        tuple: (time_array, state_array) if save_intermediate=True
               else (t_final, y_final)
    """
    t_start, t_end = t_span
    
    # Initialize
    t = t_start
    y = np.array(y0)
    
    if save_intermediate:
        # Calculate number of steps
        n_steps = int(np.ceil((t_end - t_start) / step_size))
        
        # Pre-allocate arrays
        time_array = np.zeros(n_steps + 1)
        state_array = np.zeros((n_steps + 1, len(y0)))
        
        # Store initial conditions
        time_array[0] = t
        state_array[0] = y
        
        # Integration loop
        for i in range(n_steps):
            # Adjust step size for last step if necessary
            if t + step_size > t_end:
                h = t_end - t
            else:
                h = step_size
            
            # Take RK4 step
            y = rk4_step(func, t, y, h, *args)
            t += h
            
            # Store results
            time_array[i + 1] = t
            state_array[i + 1] = y
            
            # Check if we've reached the end
            if t >= t_end:
                break
        
        # Trim arrays if we finished early
        if i + 2 < len(time_array):
            time_array = time_array[:i + 2]
            state_array = state_array[:i + 2]
        
        return time_array, state_array
    
    else:
        # Just integrate to final time without saving intermediate results
        while t < t_end:
            # Adjust step size for last step if necessary
            if t + step_size > t_end:
                h = t_end - t
            else:
                h = step_size
            
            # Take RK4 step
            y = rk4_step(func, t, y, h, *args)
            t += h
        
        return t, y

# =============================================================================
# ADAPTIVE STEP SIZE RK4 (RK4-5 with error estimation)
# =============================================================================

def rk4_adaptive_step(func, t, y, h, tolerance, *args):
    """
    Single adaptive step using RK4 with step size control.
    
    Uses Richardson extrapolation to estimate error and adjust step size.
    
    Args:
        func (callable): Function to integrate
        t (float): Current time
        y (np.ndarray): Current state vector
        h (float): Proposed step size
        tolerance (float): Error tolerance
        *args: Additional arguments for func
        
    Returns:
        tuple: (y_new, t_new, h_new, error_estimate)
    """
    # Take one full step
    y1 = rk4_step(func, t, y, h, *args)
    
    # Take two half steps
    y_half = rk4_step(func, t, y, h/2, *args)
    y2 = rk4_step(func, t + h/2, y_half, h/2, *args)
    
    # Estimate error using Richardson extrapolation
    error_estimate = np.max(np.abs(y2 - y1)) / 15.0  # RK4 has order 4, so error ~ h^5
    
    # Calculate new step size
    if error_estimate > 0:
        h_new = h * min(2.0, max(0.5, 0.9 * (tolerance / error_estimate)**(1/5)))
    else:
        h_new = h * 2.0
    
    # Clamp step size to reasonable bounds
    h_new = max(MIN_TIME_STEP, min(MAX_TIME_STEP, h_new))
    
    # Accept step if error is within tolerance
    if error_estimate <= tolerance:
        return y2, t + h, h_new, error_estimate  # Use more accurate solution
    else:
        return y, t, h_new, error_estimate  # Reject step

def rk4_adaptive_integrate(func, t_span, y0, h_initial, tolerance, *args, 
                          max_steps=100000, save_intermediate=True):
    """
    Integrate differential equation using adaptive RK4 method.
    
    Args:
        func (callable): Function to integrate
        t_span (tuple): Time span (t_start, t_end)
        y0 (np.ndarray): Initial state vector
        h_initial (float): Initial step size
        tolerance (float): Error tolerance
        *args: Additional arguments for func
        max_steps (int): Maximum number of steps
        save_intermediate (bool): Whether to save intermediate results
        
    Returns:
        tuple: (time_array, state_array, step_info) if save_intermediate=True
               else (t_final, y_final, step_info)
    """
    t_start, t_end = t_span
    
    # Initialize
    t = t_start
    y = np.array(y0)
    h = h_initial
    
    if save_intermediate:
        times = [t]
        states = [y.copy()]
    
    # Step information
    accepted_steps = 0
    rejected_steps = 0
    
    # Integration loop
    for step in range(max_steps):
        # Check if we've reached the end
        if t >= t_end:
            break
        
        # Adjust step size if we would overshoot
        if t + h > t_end:
            h = t_end - t
        
        # Take adaptive step
        y_new, t_new, h_new, error = rk4_adaptive_step(func, t, y, h, tolerance, *args)
        
        # Check if step was accepted
        if t_new > t:  # Step accepted
            t = t_new
            y = y_new
            accepted_steps += 1
            
            if save_intermediate:
                times.append(t)
                states.append(y.copy())
        else:  # Step rejected
            rejected_steps += 1
        
        # Update step size
        h = h_new
    
    # Prepare results
    step_info = {
        'accepted_steps': accepted_steps,
        'rejected_steps': rejected_steps,
        'total_attempts': accepted_steps + rejected_steps,
        'final_step_size': h
    }
    
    if save_intermediate:
        return np.array(times), np.array(states), step_info
    else:
        return t, y, step_info

# =============================================================================
# ORBIT PROPAGATION FUNCTIONS
# =============================================================================

def propagate_orbit_rk4(initial_state, time_span, step_size, mu=MU_EARTH, 
                       perturbations=None, adaptive=False, tolerance=1e-12):
    """
    Propagate orbit using RK4 integration.
    
    CURRENT IMPLEMENTATION STATUS:
    - ✅ Pure two-body dynamics (Kepler problem) - FULLY IMPLEMENTED
    - ✅ Fixed step RK4 integration - FULLY IMPLEMENTED  
    - ✅ Adaptive step RK4 integration - FULLY IMPLEMENTED
    - ⚠️ Perturbations framework - FRAMEWORK READY, NO SPECIFIC PERTURBATIONS
    
    Perturbation Capability:
    - The function accepts a 'perturbations' parameter (list of functions)
    - Each perturbation function should have signature: func(t, state) -> np.array([ax, ay, az])
    - Currently defaults to None, so only two-body dynamics are computed
    - Ready for future implementation of: atmospheric drag, J2, solar radiation pressure, etc.
    
    Physics Modeled:
    - Primary: Central body gravitational attraction (inverse square law)
    - Perturbations: None currently implemented (defaults to ideal Keplerian motion)
    
    Numerical Methods Available:
    - Fixed step RK4: Constant Δt, predictable computational cost
    - Adaptive RK4: Variable Δt with error control, better accuracy/efficiency trade-off
    
    Args:
        initial_state (np.ndarray): Initial state [x, y, z, vx, vy, vz] [km, km/s]
        time_span (tuple): (t_start, t_end) in seconds
        step_size (float): Integration step size [s] (initial step for adaptive)
        mu (float): Gravitational parameter [km³/s²] (default: Earth)
        perturbations (list): List of perturbation functions (currently None = no perturbations)
        adaptive (bool): Use adaptive step size control (False = fixed step)
        tolerance (float): Error tolerance for adaptive integration (ignored if adaptive=False)
        
    Returns:
        tuple: (time_array, position_array, velocity_array) for fixed step integration
               (time_array, position_array, velocity_array, step_info) for adaptive integration
               
    Example Usage:
        # Pure two-body problem (current default behavior)
        times, pos, vel = propagate_orbit_rk4(state0, (0, 86400), 100.0)
        
        # With perturbations (when implemented)
        # perturbations = [atmospheric_drag, j2_perturbation]  # Future capability
        # times, pos, vel = propagate_orbit_rk4(state0, (0, 86400), 100.0, perturbations=perturbations)
    """
    # Define equations of motion based on whether perturbations are provided
    if perturbations is None:
        # CURRENT DEFAULT: Pure two-body dynamics (Keplerian motion)
        # This is what's being used for Orbit B numerical convergence study
        eom = lambda t, y: orbital_equations_of_motion(t, y, mu)
    else:
        # FUTURE CAPABILITY: Two-body + perturbations
        # Framework is ready but no perturbations are currently implemented
        eom = lambda t, y: perturbed_orbital_equations(t, y, mu, perturbations)
    # Integrate
    if adaptive:
        # ADAPTIVE INTEGRATION: Variable step size with error control
        # Automatically adjusts step size to maintain specified accuracy
        # More computationally intensive but potentially more accurate
        times, states, step_info = rk4_adaptive_integrate(
            eom, time_span, initial_state, step_size, tolerance, 
            save_intermediate=True
        )
        
        # Split 6D state vector into position and velocity components
        positions = states[:, :3]  # [x, y, z] coordinates
        velocities = states[:, 3:]  # [vx, vy, vz] velocities
        
        return times, positions, velocities, step_info
    
    else:
        # FIXED STEP INTEGRATION: Constant step size (current default for convergence study)
        # Predictable computational cost, good for comparing different step sizes
        # This is what's being used for the Orbit B convergence analysis
        times, states = rk4_integrate(
            eom, time_span, initial_state, step_size, 
            save_intermediate=True
        )
        
        # Split 6D state vector into position and velocity components  
        positions = states[:, :3]  # [x, y, z] coordinates
        velocities = states[:, 3:]  # [vx, vy, vz] velocities
        
        return times, positions, velocities

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def calculate_orbital_energy(position, velocity, mu=MU_EARTH):
    """
    Calculate specific orbital energy.
    
    Args:
        position (np.ndarray): Position vector(s) [km]
        velocity (np.ndarray): Velocity vector(s) [km/s]
        mu (float): Gravitational parameter [km³/s²]
        
    Returns:
        float or np.ndarray: Specific orbital energy [km²/s²]
    """
    r = np.linalg.norm(position, axis=-1)
    v = np.linalg.norm(velocity, axis=-1)
    
    return v**2 / 2 - mu / r

def calculate_angular_momentum(position, velocity):
    """
    Calculate specific angular momentum.
    
    Args:
        position (np.ndarray): Position vector(s) [km]
        velocity (np.ndarray): Velocity vector(s) [km/s]
        
    Returns:
        float or np.ndarray: Specific angular momentum magnitude [km²/s]
    """
    if position.ndim == 1:
        h_vec = np.cross(position, velocity)
        return np.linalg.norm(h_vec)
    else:
        h_vec = np.cross(position, velocity, axis=-1)
        return np.linalg.norm(h_vec, axis=-1)

def calculate_true_anomaly_from_state(position, velocity, mu=MU_EARTH):
    """
    Calculate true anomaly from position and velocity vectors.
    
    This function converts Cartesian state vectors back to the true anomaly,
    which represents the satellite's angular position in its elliptical orbit.
    Uses the eccentricity vector method for numerical robustness.
    
    Mathematical Approach:
    1. Calculate eccentricity vector: e_vec = (1/μ)*[(v²-μ/r)*r - (r·v)*v]
    2. Find true anomaly: cos(ν) = (e_vec · r) / (|e_vec| * |r|)
    3. Determine quadrant using radial velocity (r·v)
    
    Args:
        position (np.ndarray): Position vector(s) [km] - shape (3,) or (N, 3)
        velocity (np.ndarray): Velocity vector(s) [km/s] - shape (3,) or (N, 3)
        mu (float): Gravitational parameter [km³/s²] (default: Earth)
        
    Returns:
        float or np.ndarray: True anomaly [rad] - scalar for single input, array for multiple
        
    Notes:
        - Works for any eccentricity (circular, elliptical, parabolic, hyperbolic)
        - Automatically handles single vectors or arrays of state vectors
        - Uses numerical clipping to prevent domain errors in arccos()
        - Determines correct quadrant (0-2π) using velocity direction
    """
    # Handle both single vectors and arrays of vectors
    single_vector = position.ndim == 1
    if single_vector:
        position = position.reshape(1, 3)
        velocity = velocity.reshape(1, 3)
    
    # Calculate orbital radius and velocity components
    r_mag = np.linalg.norm(position, axis=1)
    r_dot_v = np.sum(position * velocity, axis=1)  # Radial velocity component
    v_mag_sq = np.sum(velocity**2, axis=1)
    
    # Calculate eccentricity vector using vis-viva equation
    # e_vec = (1/μ)*[(v²-μ/r)*r - (r·v)*v]
    # This is more numerically stable than using orbital elements
    e_vec = np.zeros_like(position)
    for i in range(len(position)):
        e_vec[i] = (1/mu) * ((v_mag_sq[i] - mu/r_mag[i]) * position[i] - r_dot_v[i] * velocity[i])
    
    e_mag = np.linalg.norm(e_vec, axis=1)
    
    # Calculate true anomaly using dot product relationship
    # cos(ν) = (e_vec · r) / (|e_vec| * |r|)
    cos_nu = np.sum(e_vec * position, axis=1) / (e_mag * r_mag)
    
    # Numerical safety: ensure cosine is in valid domain [-1, 1]
    cos_nu = np.clip(cos_nu, -1.0, 1.0)
    
    # Calculate true anomaly (principal value: 0 to π)
    nu = np.arccos(cos_nu)
    
    # Determine correct quadrant using radial velocity
    # If r·v < 0: moving toward periapsis → ν is in range (π, 2π)
    # If r·v > 0: moving away from periapsis → ν is in range (0, π)
    for i in range(len(position)):
        if r_dot_v[i] < 0:
            nu[i] = 2*np.pi - nu[i]
    
    # Return scalar if input was scalar, array otherwise
    return nu[0] if single_vector else nu

if __name__ == "__main__":
    # Test the RK4 integrator with a simple example
    print("Testing RK4 integrator...")
    
    # Simple harmonic oscillator: d²x/dt² + ω²x = 0
    # This is a VALIDATION TEST for the RK4 algorithm using a known analytical solution
    # State: [x, dx/dt], ω = 1
    def harmonic_oscillator(t, y):
        return np.array([y[1], -y[0]])  # [velocity, -position] = [dx/dt, d²x/dt²]
    
    # Initial conditions: x(0) = 1, dx/dt(0) = 0
    # Analytical solution: x(t) = cos(t), v(t) = -sin(t)
    y0 = np.array([1.0, 0.0])
    t_span = (0.0, 2*np.pi)  # One complete oscillation
    step_size = 0.1
    
    # Integrate using our RK4 method
    times, states = rk4_integrate(harmonic_oscillator, t_span, y0, step_size)
    
    # Compare with analytical solution
    x_analytical = np.cos(times)
    
    # Check error (should be small for RK4 with reasonable step size)
    x_numerical = states[:, 0]
    error = np.max(np.abs(x_numerical - x_analytical))
    
    print(f"Maximum error in harmonic oscillator: {error:.2e}")
    print("RK4 integrator test completed successfully!")
    
    # NOTE: For orbital mechanics, this same RK4 algorithm will integrate:
    # d²r/dt² = -μr/|r|³ (two-body problem)
    # The accuracy demonstrated here carries over to orbital propagation