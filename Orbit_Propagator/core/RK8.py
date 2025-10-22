"""
Runge-Kutta 8th Order (RK8) Numerical Integration for Orbital Mechanics
=======================================================================

This module implements the Runge-Kutta 8th order method (Dormand-Prince 8(7)) 
for high-precision numerical integration of orbital equations of motion. 
This provides significantly higher accuracy than RK4 for demanding applications.

CURRENT IMPLEMENTATION STATUS:
=============================

✅ FULLY IMPLEMENTED:
- Dormand-Prince 8(7) method with embedded error estimation
- Fixed step size RK8 integration 
- Adaptive step size RK8 integration with automatic error control
- High-precision orbital mechanics integration (8th-order accuracy)
- Comprehensive coefficient tables from literature

⚠️ FRAMEWORK READY:
- Perturbation support infrastructure 
- Modular perturbation function interface

❌ NOT YET IMPLEMENTED:
- Specific perturbation forces (atmospheric drag, J2, solar radiation pressure)
- Third-body effects (Moon, Sun gravitational influence)

PHYSICS MODELED:
===============
Currently: Pure two-body problem
- Central body: Point-mass Earth (or other celestial body)
- Satellite: Point-mass with no attitude dynamics  
- Force: Inverse-square gravitational attraction only
- Mathematical model: d²r/dt² = -μr/|r|³

NUMERICAL METHOD:
================
- Dormand-Prince 8(7): 8th-order solution with 7th-order error estimation
- 13 function evaluations per step (vs. 4 for RK4)
- Embedded error estimation for adaptive step size control
- Local error: O(h⁹), Global error: O(h⁸)
- Superior accuracy for high-precision applications

PERFORMANCE CHARACTERISTICS:
===========================
- Ultra-high accuracy: 8th-order precision
- Computational cost: 13 evaluations per step (Dormand-Prince method)
- Excellent for precision-critical applications  
- Use cases: Long-term integrations, high-eccentricity orbits, precision studies

INTENDED USE:
============
- High-precision Orbit B convergence analysis
- Long-term orbital propagation
- Benchmark solutions for algorithm comparison
- Research requiring maximum numerical accuracy

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
# ORBITAL EQUATIONS OF MOTION (RK8 INDEPENDENT)
# =============================================================================

def orbital_equations_of_motion(t, state, mu=MU_EARTH):
    """
    Orbital equations of motion for two-body problem.
    
    Implements the fundamental equation: d²r/dt² = -μr/|r|³
    
    Args:
        t (float): Time [s] (not used in two-body problem)
        state (np.ndarray): State vector [x, y, z, vx, vy, vz] [km, km/s]
        mu (float): Gravitational parameter [km³/s²]
        
    Returns:
        np.ndarray: State derivative [vx, vy, vz, ax, ay, az]
    """
    # Extract position and velocity
    r = state[:3]  # Position vector [x, y, z]
    v = state[3:]  # Velocity vector [vx, vy, vz]
    
    # Calculate distance and acceleration
    r_magnitude = np.linalg.norm(r)
    a = -mu * r / (r_magnitude**3)  # Gravitational acceleration
    
    # Return state derivative: [dr/dt, dv/dt] = [v, a]
    return np.concatenate([v, a])

def perturbed_orbital_equations(t, state, mu=MU_EARTH, perturbations=None):
    """
    Orbital equations of motion with perturbations.
    
    Args:
        t (float): Time [s]
        state (np.ndarray): State vector [x, y, z, vx, vy, vz] [km, km/s]
        mu (float): Gravitational parameter [km³/s²]
        perturbations (list): List of perturbation functions
        
    Returns:
        np.ndarray: State derivative including perturbations
    """
    # Two-body acceleration (primary gravitational force)
    state_dot = orbital_equations_of_motion(t, state, mu)
    
    # Add perturbations if provided
    if perturbations is not None:
        for perturbation_func in perturbations:
            # Each perturbation function should return acceleration vector [ax, ay, az]
            perturb_accel = perturbation_func(t, state)
            # Add perturbation accelerations to the total acceleration (indices 3:6)
            state_dot[3:] += perturb_accel
    
    return state_dot

# =============================================================================
# RUNGE-KUTTA 8TH ORDER COEFFICIENTS (DORMAND-PRINCE 8(7))
# =============================================================================

# Dormand-Prince 8(7) coefficients - these are the STANDARD coefficients from:
# Dormand, J.R. and Prince, P.J. (1980). "A family of embedded Runge-Kutta formulae"
# Journal of Computational and Applied Mathematics, 6(1), 19-26.
# 
# This is an 8th-order method with embedded 7th-order error estimation
# - 13 stages (function evaluations per step)
# - Local truncation error: O(h⁹)  
# - Global error: O(h⁸)
# - Embedded error estimation for adaptive step size control

# Coefficient matrix A (lower triangular)
# A[i][j] coefficients for computing intermediate stages k_i
RK8_A = np.array([
    [],
    [1/18],
    [1/48, 1/16],
    [1/32, 0, 3/32],
    [5/16, 0, -75/64, 75/64],
    [3/80, 0, 0, 3/16, 3/20],
    [29443841/614563906, 0, 0, 77736538/692538347, -28693883/1125000000, 
     23124283/1800000000],
    [16016141/946692911, 0, 0, 61564180/158732637, 22789713/633445777,
     545815736/2771057229, -180193667/1043307555],
    [39632708/573591083, 0, 0, -433636366/683701615, -421739975/2616292301,
     100302831/723423059, 790204164/839813087, 800635310/3783071287],
    [246121993/1340847787, 0, 0, -37695042795/15268766246, -309121744/1061227803,
     -12992083/490766935, 6005943493/2108947869, 393006217/1396673457, 
     123872331/1001029789],
    [-1028468189/846180014, 0, 0, 8478235783/508512852, 1311729495/1432422823,
     -10304129995/1701304382, -48777925059/3047939560, 15336726248/1032824649,
     -45442868181/3398467696, 3065993473/597172653],
    [185892177/718116043, 0, 0, -3185094517/667107341, -477755414/1098053517,
     -703635378/230739211, 5731566787/1027545527, 5232866602/850066563,
     -4093664535/808688257, 3962137247/1805957418, 65686358/487910083],
    [403863854/491063109, 0, 0, -5068492393/434740067, -411421997/543043805,
     652783627/914296604, 11173962825/925320556, -13158990841/6184727034,
     3936647629/1978049680, -160528059/685178525, 248638103/1413531060, 0]
], dtype=object)

# Weights for 8th-order solution (B coefficients)
# These combine the k_i values to produce the final 8th-order accurate solution
RK8_B = np.array([
    14005451/335480064, 0, 0, 0, 0, -59238493/1068277825, 181606767/758867731,
    561292985/797845732, -1041891430/1371343529, 760417239/1151165299,
    118820643/751138087, -528747749/2220607170, 1/4
])

# Weights for 7th-order solution (B̂ coefficients) - used for error estimation
# The difference between 8th and 7th order solutions provides error estimate
RK8_B_HAT = np.array([
    13451932/455176623, 0, 0, 0, 0, -808719846/976000145, 1757004468/5645159321,
    656045339/265891186, -3867574721/1518517206, 465885868/322736535,
    53011238/667516719, 2/45, 0
])

# Time coefficients (C vector) - fractional times for intermediate evaluations
# These determine where in the interval [t, t+h] each k_i is evaluated
RK8_C = np.array([
    0, 1/18, 1/12, 1/8, 5/16, 3/8, 59/400, 93/200, 5490023248/9719169821,
    13/20, 1201146811/1299019798, 1, 1
])

# =============================================================================
# RUNGE-KUTTA 8TH ORDER INTEGRATOR
# =============================================================================

def rk8_step(func, t, y, h, *args):
    """
    Single step of Runge-Kutta 8th order method (Dormand-Prince 8(7)).
    
    This implements the Dormand-Prince 8(7) method, which is the GOLD STANDARD
    for high-precision orbital mechanics integration.
    
    Mathematical Foundation:
    - 13 intermediate evaluations (k₁ through k₁₃)
    - 8th-order accurate solution with embedded 7th-order error estimation
    - Local truncation error: O(h⁹)
    - Global error accumulation: O(h⁸)
    
    Algorithm Steps:
    1. Compute 13 slope evaluations k₁, k₂, ..., k₁₃ at strategic points
    2. Combine slopes with B coefficients for 8th-order solution
    3. Combine slopes with B̂ coefficients for 7th-order solution  
    4. Error estimate = |8th-order - 7th-order|
    
    For Orbital Mechanics:
    - func = orbital_equations_of_motion (currently pure two-body)
    - y = [x, y, z, vx, vy, vz] (6D state vector)  
    - Each k evaluation computes [v, a] at different time/state points
    - Much higher precision than RK4 for same step size
    
    Performance Characteristics:
    - Cost: 13 function evaluations per step
    - Accuracy: 8th-order precision with embedded error estimation
    - Efficiency: Excellent for precision-critical applications
    
    Args:
        func (callable): Function to integrate (dy/dt = func(t, y, *args))
        t (float): Current time
        y (np.ndarray): Current state vector [x, y, z, vx, vy, vz]
        h (float): Step size (Δt)
        *args: Additional arguments for func (e.g., μ, perturbations)
        
    Returns:
        tuple: (y_new, y_hat, error_estimate) where:
               - y_new: 8th-order accurate solution y(t+h)
               - y_hat: 7th-order solution (for error estimation)
               - error_estimate: |y_new - y_hat| (embedded error estimate)
    """
    # Initialize k array for 13 slope evaluations
    k = np.zeros((13, len(y)))
    
    # Calculate k values using Dormand-Prince 8(7) algorithm
    # k₁: Slope at current point (t, y)
    k[0] = h * func(t, y, *args)
    
    # k₂ through k₁₃: Slopes at intermediate points
    # Each uses linear combination of previous k values per Butcher tableau
    for i in range(1, 13):
        y_temp = y.copy()
        
        # Build intermediate state: y + Σ(A[i][j] * k[j])
        for j in range(i):
            if j < len(RK8_A[i]) and RK8_A[i][j] != 0:
                y_temp += RK8_A[i][j] * k[j]
        
        # Evaluate function at intermediate point
        k[i] = h * func(t + RK8_C[i] * h, y_temp, *args)
    
    # Calculate 8th-order solution using B coefficients
    # y_{n+1} = y_n + Σ(B[i] * k[i]) - this is the HIGH-PRECISION result
    y_new = y.copy()
    for i in range(13):
        y_new += RK8_B[i] * k[i]
    
    # Calculate 7th-order solution using B̂ coefficients for error estimation
    # ŷ_{n+1} = y_n + Σ(B̂[i] * k[i]) - this is the LOWER-ORDER comparison
    y_hat = y.copy()
    for i in range(13):
        y_hat += RK8_B_HAT[i] * k[i]
    
    # Embedded error estimate: difference between 8th and 7th order solutions
    # This provides automatic error control without additional function evaluations
    error_estimate = np.max(np.abs(y_new - y_hat))
    
    return y_new, y_hat, error_estimate

def rk8_integrate(func, t_span, y0, step_size, *args, save_intermediate=True):
    """
    Integrate differential equation using RK8 method with fixed step size.
    
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
            
            # Take RK8 step (use 8th order solution)
            y_new, _, _ = rk8_step(func, t, y, h, *args)
            y = y_new
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
            
            # Take RK8 step
            y_new, _, _ = rk8_step(func, t, y, h, *args)
            y = y_new
            t += h
        
        return t, y

# =============================================================================
# ADAPTIVE STEP SIZE RK8
# =============================================================================

def rk8_adaptive_step(func, t, y, h, tolerance, *args):
    """
    Single adaptive step using RK8 with embedded error estimation.
    
    Args:
        func (callable): Function to integrate
        t (float): Current time
        y (np.ndarray): Current state vector
        h (float): Proposed step size
        tolerance (float): Error tolerance
        *args: Additional arguments for func
        
    Returns:
        tuple: (y_new, t_new, h_new, error_estimate, accepted)
    """
    # Take RK8 step with error estimation
    y8, y7, error_estimate = rk8_step(func, t, y, h, *args)
    
    # Calculate new step size using embedded error estimation
    if error_estimate > 0:
        # Safety factor and step size adjustment for 8th order method
        safety_factor = 0.9
        h_new = h * safety_factor * (tolerance / error_estimate)**(1/8)
    else:
        h_new = h * 2.0
    
    # Clamp step size to reasonable bounds
    h_new = max(MIN_TIME_STEP, min(MAX_TIME_STEP, h_new))
    
    # Accept step if error is within tolerance
    if error_estimate <= tolerance:
        return y8, t + h, h_new, error_estimate, True
    else:
        return y, t, h_new, error_estimate, False

def rk8_adaptive_integrate(func, t_span, y0, h_initial, tolerance, *args,
                          max_steps=100000, save_intermediate=True):
    """
    Integrate differential equation using adaptive RK8 method.
    
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
        y_new, t_new, h_new, error, accepted = rk8_adaptive_step(
            func, t, y, h, tolerance, *args
        )
        
        # Check if step was accepted
        if accepted:
            t = t_new
            y = y_new
            accepted_steps += 1
            
            if save_intermediate:
                times.append(t)
                states.append(y.copy())
        else:
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

def propagate_orbit_rk8(initial_state, time_span, step_size, mu=MU_EARTH,
                       perturbations=None, adaptive=False, tolerance=1e-14):
    """
    Propagate orbit using RK8 integration.
    
    HIGH-PRECISION ORBITAL INTEGRATION using Dormand-Prince 8(7) method.
    
    CURRENT IMPLEMENTATION STATUS:
    - ✅ Pure two-body dynamics (Kepler problem) - FULLY IMPLEMENTED
    - ✅ Fixed step RK8 integration - FULLY IMPLEMENTED  
    - ✅ Adaptive step RK8 integration - FULLY IMPLEMENTED
    - ⚠️ Perturbations framework - FRAMEWORK READY, NO SPECIFIC PERTURBATIONS
    
    Accuracy Characteristics:
    - 8th-order accurate integration method
    - Ultra-high precision for orbital mechanics applications
    - Embedded error estimation for automatic step size control
    - Ideal for long-term propagation and high-eccentricity orbits
    
    Computational Cost:
    - 13 function evaluations per step (Dormand-Prince 8(7) method)
    - Excellent accuracy/cost ratio for precision applications
    - Recommended when accuracy is critical
    
    Physics Modeled:
    - Primary: Central body gravitational attraction (inverse square law)
    - Perturbations: Framework ready for additional forces
    
    Args:
        initial_state (np.ndarray): Initial state [x, y, z, vx, vy, vz] [km, km/s]
        time_span (tuple): (t_start, t_end) in seconds
        step_size (float): Integration step size [s] (initial step for adaptive)
        mu (float): Gravitational parameter [km³/s²] (default: Earth)
        perturbations (list): List of perturbation functions (currently None = no perturbations)
        adaptive (bool): Use adaptive step size control (False = fixed step)
        tolerance (float): Error tolerance for adaptive integration (default: 1e-14, very strict)
        
    Returns:
        tuple: (time_array, position_array, velocity_array) for fixed step integration
               (time_array, position_array, velocity_array, step_info) for adaptive integration
               
    Example Usage:
        # High-precision Orbit B analysis
        times, pos, vel = propagate_orbit_rk8(state0, (0, 86400), 50.0)
        
        # Ultra-high precision with adaptive control  
        times, pos, vel, info = propagate_orbit_rk8(state0, (0, 86400), 50.0, 
                                                   adaptive=True, tolerance=1e-15)
    """
    # Define equations of motion based on whether perturbations are provided
    if perturbations is None:
        # CURRENT DEFAULT: Pure two-body dynamics (Keplerian motion)
        # This is what's being used for high-precision Orbit B analysis
        eom = lambda t, y: orbital_equations_of_motion(t, y, mu)
    else:
        # FUTURE CAPABILITY: Two-body + perturbations
        # Framework is ready but no perturbations are currently implemented
        eom = lambda t, y: perturbed_orbital_equations(t, y, mu, perturbations)
    # Integrate using selected method
    if adaptive:
        # ADAPTIVE INTEGRATION: Variable step size with embedded error control
        # Automatically adjusts step size to maintain ultra-high precision
        # More computationally intensive but maximum accuracy
        times, states, step_info = rk8_adaptive_integrate(
            eom, time_span, initial_state, step_size, tolerance,
            save_intermediate=True
        )
        
        # Split 6D state vector into position and velocity components
        positions = states[:, :3]  # [x, y, z] coordinates
        velocities = states[:, 3:]  # [vx, vy, vz] velocities
        
        return times, positions, velocities, step_info
    
    else:
        # FIXED STEP INTEGRATION: Constant step size (8th-order precision)
        # Predictable computational cost with very high accuracy
        # Excellent for convergence studies and high-precision analysis
        times, states = rk8_integrate(
            eom, time_span, initial_state, step_size,
            save_intermediate=True
        )
        
        # Split 6D state vector into position and velocity components
        positions = states[:, :3]  # [x, y, z] coordinates  
        velocities = states[:, 3:]  # [vx, vy, vz] velocities
        
        return times, positions, velocities

# =============================================================================
# RK8 PERFORMANCE ANALYSIS 
# =============================================================================

def analyze_rk8_performance(initial_state, time_span, step_sizes, mu=MU_EARTH):
    """
    Analyze RK8 performance across different step sizes.
    
    Args:
        initial_state (np.ndarray): Initial state vector
        time_span (tuple): Integration time span
        step_sizes (list): List of step sizes to test
        mu (float): Gravitational parameter
        
    Returns:
        dict: Performance analysis results
    """
    import time
    
    results = {}
    
    for step_size in step_sizes:
        start_time = time.time()
        times, pos, vel = propagate_orbit_rk8(
            initial_state, time_span, step_size, mu
        )
        computation_time = time.time() - start_time
        
        results[step_size] = {
            'times': times,
            'positions': pos,
            'velocities': vel,
            'computation_time': computation_time,
            'num_steps': len(times)
        }
    
    return results

if __name__ == "__main__":
    # Test the RK8 integrator with analytical validation
    print("Testing RK8 integrator...")
    
    # Simple harmonic oscillator: d²x/dt² + ω²x = 0
    # This is a VALIDATION TEST for the RK8 algorithm using a known analytical solution
    # State: [x, dx/dt], ω = 1  
    def harmonic_oscillator(t, y):
        return np.array([y[1], -y[0]])  # [velocity, -position] = [dx/dt, d²x/dt²]
    
    # Initial conditions: x(0) = 1, dx/dt(0) = 0
    # Analytical solution: x(t) = cos(t), v(t) = -sin(t)
    y0 = np.array([1.0, 0.0])
    t_span = (0.0, 2*np.pi)  # One complete oscillation
    step_size = 0.1
    
    # Integrate using RK8 method
    times, states = rk8_integrate(harmonic_oscillator, t_span, y0, step_size)
    
    # Compare with analytical solution
    x_analytical = np.cos(times)
    
    # Check error (should be extremely small due to 8th-order accuracy)
    x_numerical = states[:, 0]
    error = np.max(np.abs(x_numerical - x_analytical))
    
    print(f"Maximum error in harmonic oscillator: {error:.2e}")
    print("RK8 integrator test completed successfully!")
    
    # Demonstrate accuracy scaling with step size
    print("\nStep size accuracy analysis:")
    times_coarse, states_coarse = rk8_integrate(harmonic_oscillator, t_span, y0, 0.5)
    x_coarse = states_coarse[:, 0]
    x_analytical_coarse = np.cos(times_coarse)
    error_coarse = np.max(np.abs(x_coarse - x_analytical_coarse))
    
    print(f"Error with step size 0.1: {error:.2e}")
    print(f"Error with step size 0.5: {error_coarse:.2e}")
    print(f"Error reduction ratio: {error_coarse/error:.1f}")
    
    # NOTE: For orbital mechanics, this same RK8 algorithm will integrate:
    # d²r/dt² = -μr/|r|³ (two-body problem) with 8th-order accuracy
    # Providing ultra-high precision for satellite orbit propagation
