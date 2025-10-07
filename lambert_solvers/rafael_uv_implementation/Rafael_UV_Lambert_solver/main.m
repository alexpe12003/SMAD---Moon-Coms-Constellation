%% main.m — Lambert (Universal Variables) Report
% Clean driver that:
%  1) defines (r1, r2, dt, mu)
%  2) runs lambert_uv
%  3) computes COEs from (r1, v1)
%  4) prints a tidy report
%
% Notes:
%  • Propagation (kepler_universal) and plotting are intentionally omitted.
% 

clear; clc;

%% 1) Problem definition (consistent units)
mu = 398600.4418;                % [km^3/s^2] Earth
r1 = [  5000; 10000; 2100 ];     % [km]
r2 = [ -14600;  2500; 7000 ];    % [km]
dt = 1*60*60;                    % [s] 1 hour

%% 2) Run the Universal-Variables Lambert solver
% Adjust options as needed:
%   'z0',1.5         -> seed for z (your preferred initial guess)
%   'longway',true   -> use long-way arc (Δθ in (π, 2π))
%   'tol',1e-10      -> residual tolerance |F(z)| in seconds
%   'maxIter',50     -> iteration cap
%   'verbose',true   -> iteration trace
[v1, v2, info] = lambert_uv(r1, r2, dt, mu, 'z0', 1.5, 'verbose', true);

%% 3) Compute classical orbital elements from (r1, v1)
% coe = [h, e, RAAN, i, argPeri, trueAnom, a]
coe = coe_from_sv(r1, v1, mu);

%% 4) Report
fprintf('\n================ LAMBERT (UV) REPORT ================\n');
fprintf('r1 = [%10.3f %10.3f %10.3f] km\n', r1);
fprintf('r2 = [%10.3f %10.3f %10.3f] km\n', r2);
fprintf('Δt = %g s (%.3f h)   μ = %.6f km^3/s^2\n\n', dt, dt/3600, mu);

% Convergence block
fprintf('Method: Universal Variables (UV)\n');
fprintf('  converged: %s   iters: %d   |F|: %.3e s   z: %.6g\n', ...
    yesno(getfielddef(info,'converged',false)), ...
    getfielddef(info,'iterations',NaN), ...
    getfielddef(info,'tof_err',NaN), ...
    getfielddef(info,'z',NaN));
if isfield(info,'message')
    fprintf('  note: %s\n', info.message);
end

% Velocities
fprintf('  v1 (km/s): [%12.6f %12.6f %12.6f]\n', v1);
fprintf('  v2 (km/s): [%12.6f %12.6f %12.6f]\n\n', v2);

% COEs (angles in degrees)
h     = coe(1);
e     = coe(2);
RAAN  = rad2deg(coe(3));
inc   = rad2deg(coe(4));
argp  = rad2deg(coe(5));
nu    = rad2deg(coe(6));
a     = coe(7);

fprintf('  COE (from r1,v1):\n');
fprintf('    h      = % .6f km^2/s\n', h);
fprintf('    e      = % .8f [-]\n',   e);
fprintf('    RAAN   = % .6f deg\n',   RAAN);
fprintf('    i      = % .6f deg\n',   inc);
fprintf('    ω      = % .6f deg\n',   argp);
fprintf('    ν      = % .6f deg\n',   nu);
fprintf('    a      = % .6f km\n',    a);
fprintf('=====================================================\n');

%% LATER: Propagation check (endpoint error)
% Add kepler_universal.m, you can verify the solution:
% [r2p, v2p] = kepler_universal(r1, v1, dt, mu);
% fprintf('Endpoint error |r2p - r2| = %.6e km\n', norm(r2p - r2));

%% LATER: Plotting
% After you add plot_transfer3d.m, you can visualize the arc:
% plot_transfer3d(r1, r2, v1, mu, dt, 'BodyRadius', 6378, 'BodyName', 'Earth');

%% ----------------- helpers (local) -----------------
function s = yesno(b)
% Map logical to 'yes'/'no' for neat printing
    if islogical(b) && isscalar(b)
        s = ternary(b,'yes','no');
    else
        s = 'n/a';
    end
end

function y = ternary(cond, a, b)
    if cond, y = a; else, y = b; end
end

function v = getfielddef(S, fld, default)
% get struct field with a default if missing
    if isstruct(S) && isfield(S, fld), v = S.(fld); else, v = default; end
end

%% Sanity checks