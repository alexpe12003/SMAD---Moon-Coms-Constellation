function [R, V, info] = kepler_universal(R0, V0, dt, mu, varargin)
% KEPLER_UNIVERSAL  Two-body propagation using universal variables.
%
%   [R, V, info] = kepler_universal(R0, V0, dt, mu, 'Name',Value,...)
%
% Inputs
%   R0  : 3x1 initial position [km]
%   V0  : 3x1 initial velocity [km/s]
%   dt  : propagation time     [s]  (can be positive or negative)
%   mu  : gravitational parameter [km^3/s^2]
%
% Name–Value options
%   'tol'     (1e-10) : convergence on |F(chi)| (seconds)
%   'maxIter' (50)    : max Newton iterations
%   'verbose' (false) : print iteration trace
%
% Outputs
%   R, V : 3x1 position and velocity at t = dt
%   info : diagnostics (converged, iterations, chi, z, message)
%
% Notes
%   • Works for elliptic, parabolic (α→0), and hyperbolic cases.
%   • Requires STUMPFF.M providing C(z), S(z).
%   • Reference equations follow standard Vallado universal formulation.

% ---------- parse options ----------
ip = inputParser;
addParameter(ip,'tol',1e-10,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip,'maxIter',50,@(x)isnumeric(x)&&isscalar(x)&&x>=5);
addParameter(ip,'verbose',false,@(x)islogical(x)||ismember(x,[0,1]));
parse(ip,varargin{:});
opt = ip.Results;

% ---------- input hygiene ----------
R0 = R0(:); V0 = V0(:);
if numel(R0)~=3 || numel(V0)~=3, error('R0 and V0 must be 3x1 vectors.'); end
if ~isscalar(mu) || mu<=0, error('mu must be a positive scalar.'); end
if ~isscalar(dt), error('dt must be a scalar.'); end

% ---------- scalars from initial state ----------
r0  = norm(R0);
v0  = norm(V0);
vr0 = dot(R0,V0)/r0;

% alpha = 2/r0 - v0^2/mu  (reciprocal of a with sign)
alpha = 2/r0 - (v0^2)/mu;      % [1/km] ; alpha>0 elliptic, =0 parabola, <0 hyperbolic

% ---------- choose a starter for universal anomaly chi ----------
% Units of chi are [sqrt(mu)*s] i.e. km^(3/2)/s
if alpha > 0
    chi = sqrt(mu)*alpha*dt;                  % good starter for elliptic
elseif alpha < 0
    % hyperbolic: a simple robust starter proportional to |dt|
    chi = sign(dt) * sqrt(-1/alpha) * log( (-2*mu*alpha*abs(dt)) / ...
         (vr0 + sign(dt)*sqrt(-mu*alpha)*(1 - r0*alpha)) );
    if ~isfinite(chi)
        chi = sign(dt) * sqrt(mu) * abs(alpha) * dt;  % fallback
    end
else
    % near parabolic
    chi = sqrt(mu) * dt / r0;                 % Ok for α≈0
end

% Guard against NaN/Inf starters
if ~isfinite(chi), chi = sqrt(mu)*alpha*dt; end
if ~isfinite(chi), chi = 0;                 end

% ---------- Newton iterations on F(chi) = 0 ----------
% F(chi) = (r0*vr0/sqrt(mu)) * chi^2 * C(z) + (1 - alpha*r0) * chi^3 * S(z) + r0*chi - sqrt(mu)*dt
% with z = alpha * chi^2
iter = 0; converged = false; msg = 'OK';

while iter < opt.maxIter
    iter = iter + 1;

    z = alpha * chi^2;
    [C, S] = stumpff(z);

    % F(chi)
    F = (r0*vr0/sqrt(mu))*chi^2*C + (1 - alpha*r0)*chi^3*S + r0*chi - sqrt(mu)*dt;

    if opt.verbose
        fprintf('it=%2d  chi=% .6e  F=% .3e  z=% .3e\n', iter, chi, F, z);
    end

    if abs(F) < opt.tol
        converged = true; break;
    end

    % dF/dchi (Vallado, stable form)
    dF = (r0*vr0/sqrt(mu))*(1 - z*S) + (1 - alpha*r0)*chi*(1 - z*C) + r0;

    % Newton step with a little safeguarding
    if ~isfinite(dF) || dF == 0
        step = sign(F) * max(1, abs(F));   % crude fallback
    else
        step = F / dF;
    end

    chi_new = chi - step;

    % Gentle damping if the step is too large (helps in nasty cases)
    if ~isfinite(chi_new)
        chi_new = 0.5*chi;
    elseif abs(chi_new) > 10*max(1,abs(chi))
        chi_new = chi - 0.25*step;
    end

    chi = chi_new;
end

if ~converged
    msg = 'Did not converge';
end

% ---------- Lagrange f, g, fdot, gdot and final state ----------
z    = alpha * chi^2;
[C,S] = stumpff(z);

f    = 1 - (chi^2/r0)*C;
g    = dt - (1/sqrt(mu))*chi^3*S;

R = f*R0 + g*V0;
r = norm(R);

% fdot = sqrt(mu)/(r*r0) * (alpha*chi^3*S - chi)
% gdot = 1 - (chi^2/r) * C
fdot = (sqrt(mu)/(r*r0)) * (alpha*chi^3*S - chi);
gdot = 1 - (chi^2/r) * C;

V = fdot*R0 + gdot*V0;

% ---------- package info ----------
info = struct('converged',converged, 'iterations',iter, ...
              'chi',chi, 'z',z, 'message',msg);
end
