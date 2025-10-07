function coe = coe_from_sv(R, V, mu)
% COE_FROM_SV  Classical Orbital Elements from state vectors.
%
%   coe = coe_from_sv(R, V, mu)
%
% Inputs
%   R  : 3x1 position vector in an inertial, Earth-centered frame [km]
%   V  : 3x1 velocity vector in the same frame                    [km/s]
%   mu : gravitational parameter                                   [km^3/s^2]
%
% Output
%   coe : 1x7 vector of classical orbital elements (radians)
%         [ h, e, RAAN, i, argPeri, trueAnom, a ]
%           h        : specific angular momentum magnitude [km^2/s]
%           e        : scalar eccentricity [-]
%           RAAN     : right ascension of ascending node, Ω [rad]
%           i        : inclination, i [rad]
%           argPeri  : argument of perigee, ω [rad]
%           trueAnom : true anomaly, ν [rad]
%           a        : semi-major axis [km]  (a<0 for hyperbola)
%
% Notes
%   • Angles are in radians.
%   • Robust to circular & equatorial special cases: Ω and ω are set to 0
%     when undefined (n≈0 or e≈0).
%   • Uses numerically safe acos by clamping inputs to [-1,1].

    % ---------- Input hygiene ----------
    R = R(:);  V = V(:);
    if numel(R) ~= 3 || numel(V) ~= 3
        error('coe_from_sv:badInput','R and V must be 3-element vectors.');
    end
    if ~isscalar(mu) || mu <= 0
        error('coe_from_sv:badMu','mu must be a positive scalar.');
    end

    % Small thresholds (tune if needed)
    eps_e = 1e-10;   % eccentricity ~ 0 threshold
    eps_n = 1e-12;   % node vector magnitude ~ 0 threshold

    % Utility for safe acos
    clamp = @(x) max(-1,min(1,x));

    % ---------- Scalars from state ----------
    r  = norm(R);
    v  = norm(V);
    vr = dot(R,V) / r;                    % radial velocity component [km/s]

    % ---------- Specific angular momentum ----------
    H = cross(R,V);                        % vector [km^2/s]
    h = norm(H);

    % ---------- Inclination (i) ----------
    % Guard acos domain with clamp to handle roundoff
    i = acos( clamp( H(3) / h ) );

    % ---------- Node vector (line of nodes) ----------
    K = [0;0;1];                           % reference z-axis
    N = cross(K, H);                       % points toward ascending node
    n = norm(N);

    % ---------- RAAN (Ω) ----------
    if n > eps_n
        % Ω = acos(N_x / |N|); if N_y < 0, take 2π − Ω
        RAAN = acos( clamp( N(1) / n ) );
        if N(2) < 0
            RAAN = 2*pi - RAAN;
        end
    else
        % Equatorial orbit: node line undefined
        RAAN = 0;
    end

    % ---------- Eccentricity vector & magnitude ----------
    % e_vec = [ (v^2 - mu/r) R - (r*vr) V ] / mu
    E_vec = ( (v^2 - mu/r) * R - (r*vr) * V ) / mu;
    e = norm(E_vec);

    % ---------- Argument of perigee (ω) ----------
    if n > eps_n && e > eps_e
        % ω = acos( (N·e) / (|N| e) ); if e_z < 0, take 2π − ω
        argPeri = acos( clamp( dot(N,E_vec)/(n*e) ) );
        if E_vec(3) < 0
            argPeri = 2*pi - argPeri;
        end
    else
        % Circular or equatorial: ω undefined → set to 0
        argPeri = 0;
    end

    % ---------- True anomaly (ν) ----------
    if e > eps_e
        % ν = acos( (e·R) / (e r) ); if vr < 0, take 2π − ν
        trueAnom = acos( clamp( dot(E_vec,R)/(e*r) ) );
        if vr < 0
            trueAnom = 2*pi - trueAnom;
        end
    else
        % Circular case: use node line to define ν
        if n > eps_n
            % If cp = N×R has positive z, keep acos; else wrap 2π − acos
            cp = cross(N, R);
            base = acos( clamp( dot(N,R)/(n*r) ) );
            if cp(3) >= 0
                trueAnom = base;
            else
                trueAnom = 2*pi - base;
            end
        else
            % Circular & equatorial: ν undefined → set to 0
            trueAnom = 0;
        end
    end

    % ---------- Semi-major axis (a) ----------
    % a = h^2 / mu / (1 − e^2)   (negative for hyperbola, e>1)
    a = (h^2) / mu / (1 - e^2);

    % ---------- Pack result ----------
    coe = [h, e, RAAN, i, argPeri, trueAnom, a];
end