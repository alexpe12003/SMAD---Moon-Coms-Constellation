import numpy as np

def kepler_universal(R0, V0, dt, mu, tol=1e-10, maxIter=50, verbose=False):
    """
    Two-body propagation using universal variables.

    Parameters
    ----------
    R0 : array_like
        Initial position vector [km]
    V0 : array_like
        Initial velocity vector [km/s]
    dt : float
        Propagation time [s]
    mu : float
        Gravitational parameter [km^3/s^2]
    tol : float, optional
        Convergence tolerance on |F(chi)| (seconds)
    maxIter : int, optional
        Maximum Newton iterations
    verbose : bool, optional
        Print iteration trace

    Returns
    -------
    R : ndarray
        Position at t = dt [km]
    V : ndarray
        Velocity at t = dt [km/s]
    info : dict
        Diagnostics (converged, iterations, chi, z, message)
    """
    R0 = np.asarray(R0).reshape(3)
    V0 = np.asarray(V0).reshape(3)
    if R0.size != 3 or V0.size != 3:
        raise ValueError("R0 and V0 must be 3-element vectors.")
    if not np.isscalar(mu) or mu <= 0:
        raise ValueError("mu must be a positive scalar.")
    if not np.isscalar(dt):
        raise ValueError("dt must be a scalar.")

    r0 = np.linalg.norm(R0)
    v0 = np.linalg.norm(V0)
    vr0 = np.dot(R0, V0) / r0

    alpha = 2 / r0 - (v0 ** 2) / mu

    # Starter for universal anomaly chi
    if alpha > 0:
        chi = np.sqrt(mu) * alpha * dt
    elif alpha < 0:
        try:
            chi = np.sign(dt) * np.sqrt(-1 / alpha) * np.log(
                (-2 * mu * alpha * abs(dt)) /
                (vr0 + np.sign(dt) * np.sqrt(-mu * alpha) * (1 - r0 * alpha))
            )
        except Exception:
            chi = np.sign(dt) * np.sqrt(mu) * abs(alpha) * dt
        if not np.isfinite(chi):
            chi = np.sign(dt) * np.sqrt(mu) * abs(alpha) * dt
    else:
        chi = np.sqrt(mu) * dt / r0

    if not np.isfinite(chi):
        chi = np.sqrt(mu) * alpha * dt
    if not np.isfinite(chi):
        chi = 0

    iter = 0
    converged = False
    msg = 'OK'

    while iter < maxIter:
        iter += 1
        z = alpha * chi ** 2
        C, S = stumpff(z)  # You must implement stumpff(z) to return C(z), S(z)

        F = (r0 * vr0 / np.sqrt(mu)) * chi ** 2 * C + \
            (1 - alpha * r0) * chi ** 3 * S + r0 * chi - np.sqrt(mu) * dt

        if verbose:
            print(f"it={iter:2d}  chi={chi: .6e}  F={F: .3e}  z={z: .3e}")

        if abs(F) < tol:
            converged = True
            break

        dF = (r0 * vr0 / np.sqrt(mu)) * (1 - z * S) + \
             (1 - alpha * r0) * chi * (1 - z * C) + r0

        if not np.isfinite(dF) or dF == 0:
            step = np.sign(F) * max(1, abs(F))
        else:
            step = F / dF

        chi_new = chi - step

        if not np.isfinite(chi_new):
            chi_new = 0.5 * chi
        elif abs(chi_new) > 10 * max(1, abs(chi)):
            chi_new = chi - 0.25 * step

        chi = chi_new

    if not converged:
        msg = 'Did not converge'

    z = alpha * chi ** 2
    C, S = stumpff(z)

    f = 1 - (chi ** 2 / r0) * C
    g = dt - (1 / np.sqrt(mu)) * chi ** 3 * S

    R = f * R0 + g * V0
    r = np.linalg.norm(R)

    fdot = (np.sqrt(mu) / (r * r0)) * (alpha * chi ** 3 * S - chi)
    gdot = 1 - (chi ** 2 / r) * C

    V = fdot * R0 + gdot * V0

    info = {
        'converged': converged,
        'iterations': iter,
        'chi': chi,
        'z': z,
        'message': msg
    }

    return R, V, info