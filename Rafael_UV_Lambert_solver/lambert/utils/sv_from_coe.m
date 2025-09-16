function [r, v] = sv_from_coe(coe, mu)
% SV_FROM_COE  State vector (r,v) from classical orbital elements.
%
%   [r, v] = sv_from_coe(coe, mu)
%
% Inputs
%   coe : [h e RAAN i w TA]
%         h    : specific angular momentum [km^2/s]
%         e    : eccentricity [-]
%         RAAN : right ascension of ascending node Ω [rad]
%         i    : inclination [rad]
%         w    : argument of perigee [rad]
%         TA   : true anomaly [rad]
%   mu  : gravitational parameter [km^3/s^2]
%
% Outputs
%   r : 1x3 position vector in inertial (ECI) frame [km]
%   v : 1x3 velocity vector in inertial (ECI) frame [km/s]
%
% Notes
%   - Implements Algorithm 4.2 (Curtis, Vallado).
%   - Perifocal frame (PQW) → geocentric equatorial transformation.
%   - Angles must be in radians.
%
% User M-functions required: none

% ---------- unpack orbital elements ----------
h     = coe(1);
e     = coe(2);
RAAN  = coe(3);
incl  = coe(4);
w     = coe(5);
TA    = coe(6);

% ---------- position & velocity in perifocal frame ----------
rp = (h^2/mu) * (1/(1 + e*cos(TA))) * ...
     ( cos(TA)*[1;0;0] + sin(TA)*[0;1;0] );

vp = (mu/h) * ( -sin(TA)*[1;0;0] + (e + cos(TA))*[0;1;0] );

% ---------- rotation matrices ----------
% R3(-RAAN) * R1(-incl) * R3(-w) = Q_pX
R3_W = [ cos(RAAN)  sin(RAAN) 0;
        -sin(RAAN)  cos(RAAN) 0;
              0           0   1];

R1_i = [1      0          0;
        0  cos(incl)  sin(incl);
        0 -sin(incl)  cos(incl)];

R3_w = [ cos(w)  sin(w) 0;
        -sin(w)  cos(w) 0;
           0       0    1];

Q_pX = R3_W' * R1_i' * R3_w';   % PQW → inertial

% ---------- transform to inertial frame ----------
r = (Q_pX * rp)';   % row vector
v = (Q_pX * vp)';   % row vector

end
