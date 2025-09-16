function [v1, v2, info] = lambert_uv(r1, r2, dt, mu, varargin)
%LAMBERT_UV  Universal-Variables iterative Lambert solver (0-rev).
%   Solves the classical Lambert boundary-value problem using the
%   universal-variables formulation. Iterates on the scalar z such that
%   T(z) = dt, then recovers velocities via Lagrange coefficients.
%
%   [v1, v2, info] = lambert_uv(r1, r2, dt, mu, 'Name', Value, ...)
%
% Inputs (SI-like but consistent units)
%   r1, r2 : 3x1 position vectors [km]
%   dt     : time of flight        [s]  (dt > 0)
%   mu     : gravitational parameter [km^3/s^2]
%
% Name–Value options (all optional)
%   'longway'  (false)  : if true, use the long-way transfer (Δθ ∈ (π, 2π))
%   'tol'      (1e-10)  : convergence on |F(z)| = |T(z)/sqrt(mu) - dt| [s]
%   'maxIter'  (50)     : maximum Newton steps
%   'verbose'  (false)  : print iteration trace (it, z, F)
%   'zMin'     (-4)     : lower bound for z search (hyperbolic side)
%   'zMax'     (+4)     : upper bound for z search (elliptic side)
%   'z0'       (0)      : initial guess for z (e.g., 1.5 if you expect the root there)
%
% Outputs
%   v1, v2 : 3x1 velocities at r1 and r2 [km/s]; NaNs if not converged
%   info   : struct with diagnostics
%       .converged, .iterations, .z, .tof_err
%       .theta, .A, .y, .f, .g, .gdot
%       .message
%
% Notes
%   • This routine targets 0-revolution solutions (no full wraps).
%   • Requires STUMPFF.M (C(z), S(z)) on your path.
%   • References: Vallado, Lancaster & Blanchard (UV Lambert); the dF/dz
%     formula includes the special-case z=0 expression you provided.

%----------------------------------------------------------------------
% 0) Parse & validate inputs (inputParser works on older MATLABs)
%----------------------------------------------------------------------

ip = inputParser;
ip.FunctionName = 'lambert_uv';

addRequired(ip, 'r1', @(x)isnumeric(x)&&isvector(x)&&numel(x)==3);
addRequired(ip, 'r2', @(x)isnumeric(x)&&isvector(x)&&numel(x)==3);
addRequired(ip, 'dt', @(x)isnumeric(x)&&isscalar(x)&&x>0);
addRequired(ip, 'mu', @(x)isnumeric(x)&&isscalar(x)&&x>0);

addParameter(ip, 'longway', false,  @(x)islogical(x)||ismember(x,[0,1]));
addParameter(ip, 'tol',     1e-10,  @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip, 'maxIter', 50,     @(x)isnumeric(x)&&isscalar(x)&&x>=5);
addParameter(ip, 'verbose', false,  @(x)islogical(x)||ismember(x,[0,1]));
addParameter(ip, 'zMin',    -4,     @(x)isnumeric(x)&&isscalar(x)&&x<0);
addParameter(ip, 'zMax',    +4,     @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip, 'z0',       0,     @(x)isnumeric(x)&&isscalar(x));

parse(ip, r1, r2, dt, mu, varargin{:});
opt = ip.Results;

% Normalize vector shapes and magnitudes
r1  = r1(:);   r2  = r2(:);
r1n = norm(r1); r2n = norm(r2);

if r1n==0 || r2n==0,   error('lambert_uv:badInput','r1 and r2 must be nonzero.'); end
if opt.zMin >= opt.zMax, error('lambert_uv:badBounds','zMin must be < zMax.');     end

%----------------------------------------------------------------------
% 1) Geometry: transfer angle θ and A parameter
%----------------------------------------------------------------------

% Unsigned angle via atan2 (robust to roundoff); (0, π]
cr     = cross(r1, r2);
sin_d  = norm(cr) / (r1n * r2n);                            % ≥ 0
cos_d  = max(-1, min(1, dot(r1, r2)/(r1n*r2n)));            % clamp for safety
theta  = atan2(sin_d, cos_d);                               % (0, π]

% Long-way: reflect to (π, 2π)
if opt.longway
    theta = 2*pi - theta;                                   % (π, 2π)
end

% A encodes chord/semiperimeter geometry and the short/long arc via sin(θ)
den = 1 - cos(theta);
if den <= eps
    % θ≈0 : degenerate geometry for this formulation
    v1 = nan(3,1); v2 = nan(3,1);
    info = packInfo(false, 0, NaN, Inf, theta, NaN, NaN, NaN, NaN, NaN, 'Δθ≈0: degenerate geometry');
    return;
end

A = sin(theta) * sqrt( (r1n * r2n) / den );
if abs(A) < 1e-14
    v1 = nan(3,1); v2 = nan(3,1);
    info = packInfo(false, 0, NaN, Inf, theta, A, NaN, NaN, NaN, NaN, 'A≈0: degenerate geometry');
    return;
end

%----------------------------------------------------------------------
% 2) Residual F(z) = T(z)/sqrt(mu) - dt   (nested function)
%    T(z) is UV time of flight; y(z) is the UV auxiliary variable
%----------------------------------------------------------------------

    function [F, y, C, S] = tof_residual(z)
        % Stumpff functions; series-safe near z≈0 in your stumpff.m
        [C, S] = stumpff(z);

        % Guard: for z>0 (elliptic), C should be >0; keep it positive
        if (z > 0) && (C <= 0)
            C = max(C, realmin);
        end

        % y(z) (Lancaster–Blanchard / Vallado)
        sqrtC = sqrt(C);
        y = r1n + r2n + A * (z*S - 1) / sqrtC;

        % Infeasible branch if y≤0 (no physical conic)
        if ~(isfinite(y)) || y <= 0
            F = +Inf;    % signal invalid
            return;
        end

        % Time of flight (scaled), then subtract your requested dt
        chi = sqrt(y / C);               % universal anomaly
        Tz  = chi^3 * S + A * sqrt(y);   % has units sqrt(mu)*seconds
        F   = Tz / sqrt(mu) - dt;        % residual in seconds
    end

%----------------------------------------------------------------------
% 3) Analytic derivative dF/dz (nested function)
%    Uses your provided formulas; special case at z = 0.
%----------------------------------------------------------------------

    function dF = dFdz_analytic(z)
        % Evaluate at current z
        [~, y, C, S] = tof_residual(z);
        if ~isfinite(y) || y <= 0
            dF = NaN; return;
        end

        % Switch to z=0 formula near the singular term 1/(2z)
        z_eps = 1e-12;

        if abs(z) < z_eps
            % Special case at z=0 (use y evaluated exactly at 0)
            [~, y0] = tof_residual(0.0);
            if ~isfinite(y0) || y0 <= 0
                dF = NaN; return;
            end

            % dT/dz(0) from the formula you pasted
            dTdz = (sqrt(2)/40) * (y0^(3/2)) + (A/8) * ( sqrt(y0) + A * (1/(2*y0)) );
            dF   = dTdz / sqrt(mu);
            return;
        end

        % General case z ≠ 0
        % dT/dz =
        %   (y/C)^(3/2) * { (1/(2z)) [ C - (3/2)(S/C) ] + (3/4)(S^2/C) }
        % + (A/8) * [ 3(S/C)√y + A√(C/y) ]
        % then dF/dz = dT/dz / √mu

        if C <= 0
            C = max(C, realmin);
        end

        yc        = y / C;
        sqrt_yc   = sqrt(yc);
        sqrt_y    = sqrt(y);

        curly     = (1/(2*z)) * ( C - (3/2)*(S/C) ) + (3/4)*(S*S)/C;
        term1     = (sqrt_yc^3) * curly;
        term2     = (A/8) * ( 3*(S/C)*sqrt_y + A*sqrt(C/y) );

        dTdz      = term1 + term2;
        dF        = dTdz / sqrt(mu);
    end

%----------------------------------------------------------------------
% 4) Initial guess for z (seed -> local probe -> coarse sweep)
%----------------------------------------------------------------------

z = opt.z0;                        % user-provided seed (e.g., 1.5)
[F, y] = tof_residual(z);

if ~isfinite(F)
    % (a) Small local wiggle around z0 to find a valid F
    deltas = [1e-3, 1e-2, 1e-1, 0.3, 0.7, 1.0];
    signs  = [+1, -1];
    got = false;

    for d = deltas
        for sgn = signs
            z_try = z + sgn*d;
            if z_try < opt.zMin || z_try > opt.zMax, continue; end
            [F_try, y_try] = tof_residual(z_try);
            if isfinite(F_try)
                z = z_try; F = F_try; y = y_try; got = true; break;
            end
        end
        if got, break; end
    end

    % (b) If local probe failed, coarse probe across [zMin, zMax]
    if ~got
        probes = linspace(opt.zMin, opt.zMax, nine_points(opt.zMin, opt.zMax));
        for p = probes
            [Fp, yp] = tof_residual(p);
            if isfinite(Fp)
                z = p; F = Fp; y = yp; got = true; break;
            end
        end
    end

    % (c) If still nothing valid, abort early
    if ~got
        v1 = nan(3,1); v2 = nan(3,1);
        info = packInfo(false, 0, NaN, Inf, theta, A, NaN, NaN, NaN, NaN, 'Failed to find valid initial z');
        return;
    end
end

% helper: choose 9 points unless the interval is tiny
    function n = nine_points(a, b)
        if b - a < 2
            n = 11;
        else
            n = 9;
        end
    end

%----------------------------------------------------------------------
% 5) Iterate on z (Newton with analytic derivative, safeguarded)
%----------------------------------------------------------------------

iter       = 0;
converged  = false;
last_z     = z;
last_F     = F;

while iter < opt.maxIter
    iter = iter + 1;

    % Analytic derivative dF/dz (with z=0 special case)
    dF = dFdz_analytic(z);

    % Fallback if derivative unusable (rare, but be safe)
    if ~isfinite(dF) || dF == 0
        if iter > 1 && (z ~= last_z) && isfinite(last_F)
            dF = (F - last_F) / (z - last_z);   % secant slope
        else
            dF = sign(F) * max(1, abs(F));      % crude but safe
        end
    end

    % Newton step
    z_new = z - F/dF;

    % Safeguards: keep inside [zMin,zMax] and avoid NaN/Inf jumps
    if z_new < opt.zMin, z_new = 0.5*(z + opt.zMin); end
    if z_new > opt.zMax, z_new = 0.5*(z + opt.zMax); end
    if ~isfinite(z_new), z_new = 0.5*(opt.zMin + opt.zMax); end

    % Advance and evaluate
    last_z = z; last_F = F;
    z = z_new;
    [F, y] = tof_residual(z);

    if opt.verbose
        fprintf('it=%2d  z=% .6e  F=% .3e\n', iter, z, F);
    end

    % Converged on residual (seconds)
    if isfinite(F) && abs(F) < opt.tol
        converged = true;
        break;
    end
end

% If not converged or invalid y, exit gracefully
if ~converged || ~isfinite(F) || y <= 0
    v1 = nan(3,1); v2 = nan(3,1);
    info = packInfo(false, iter, z, abs(F), theta, A, y, NaN, NaN, NaN, 'Did not converge or invalid y');
    return;
end

%----------------------------------------------------------------------
% 6) Recover velocities via Lagrange f, g, gdot
%----------------------------------------------------------------------

f    = 1 - y / r1n;
g    = A * sqrt(y / mu);
gdot = 1 - y / r2n;

if abs(g) < 1e-14 || ~isfinite(g)
    v1 = nan(3,1); v2 = nan(3,1);
    info = packInfo(false, iter, z, abs(F), theta, A, y, f, g, gdot, 'g≈0 (singular)');
    return;
end

v1 = (r2 - f*r1)   / g;
v2 = (gdot*r2 - r1)/ g;

info = packInfo(true, iter, z, abs(F), theta, A, y, f, g, gdot, 'OK');

end % <<< end main function

%======================================================================
% 7) Small helper to package diagnostics (separate/local function)
%======================================================================
function s = packInfo(conv, it, z, toferr, theta, A, y, f, g, gdot, msg)
s = struct( ...
    'converged', logical(conv), ...
    'iterations', it, ...
    'z',         z, ...
    'tof_err',   toferr, ...
    'theta',     theta, ...
    'A',         A, ...
    'y',         y, ...
    'f',         f, ...
    'g',         g, ...
    'gdot',      gdot, ...
    'message',   msg);
end


