function plot_transfer3d(r1, r2, v1, mu, dt, varargin)
% PLOT_TRANSFER3D  3D visualization of a Lambert transfer.
%
%   plot_transfer3d(r1, r2, v1, mu, dt)
%   plot_transfer3d(..., 'N',400, 'BodyRadius',6378, 'BodyName','Earth', ...
%                           'ShowVel',true, 'ShowPlane',true, 'NewFigure',true)
%
% Inputs
%   r1, r2 : 3x1 vectors, initial/target position [km]
%   v1     : 3x1 vector, initial velocity from Lambert [km/s]
%   mu     : gravitational parameter [km^3/s^2]
%   dt     : time of flight [s]
%
% Options (Name–Value)
%   'N'          (300)   number of points along the arc
%   'BodyRadius' ([])    central body radius (km) to draw a sphere
%   'BodyName'   ('')    label for the central body
%   'ShowVel'    (true)  draw v1 and terminal velocity arrows
%   'ShowPlane'  (false) draw a translucent disk of the transfer plane
%   'NewFigure'  (true)  open a new figure
%
% Notes
%   • Requires kepler_universal.m and stumpff.m on path.
%   • All units are km (positions) and km/s (velocities).

% --------- options ----------
ip = inputParser;
addParameter(ip,'N',300,@(x)isnumeric(x)&&isscalar(x)&&x>=20);
addParameter(ip,'BodyRadius',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>0));
addParameter(ip,'BodyName','',@(x)ischar(x)||isstring(x));
addParameter(ip,'ShowVel',true,@(x)islogical(x)||ismember(x,[0,1]));
addParameter(ip,'ShowPlane',false,@(x)islogical(x)||ismember(x,[0,1]));
addParameter(ip,'NewFigure',true,@(x)islogical(x)||ismember(x,[0,1]));
parse(ip,varargin{:});
opt = ip.Results;

r1 = r1(:); r2 = r2(:); v1 = v1(:);

% --------- sample times & propagate arc ----------
tvec = linspace(0, dt, opt.N);
R = zeros(3, opt.N);
V = zeros(3, opt.N);
for k = 1:opt.N
    [R(:,k), V(:,k)] = kepler_universal(r1, v1, tvec(k), mu);
end

% --------- figure setup ----------
if opt.NewFigure, figure('Color','w'); end
hold on; grid on; axis equal;
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');
title('3D Lambert Transfer');

% Axes triad
Ltri = 0.15 * max(range(R,2)); if ~isfinite(Ltri)||Ltri<=0, Ltri = 5000; end
quiver3(0,0,0, Ltri,0,0, 0,'k','LineWidth',1,'DisplayName','X');
quiver3(0,0,0, 0,Ltri,0, 0,'k','LineWidth',1,'HandleVisibility','off');
quiver3(0,0,0, 0,0,Ltri, 0,'k','LineWidth',1,'HandleVisibility','off');

% Central body (optional)
if ~isempty(opt.BodyRadius)
    draw_body(opt.BodyRadius, char(opt.BodyName));
end

% Transfer plane disk (optional)
if opt.ShowPlane
    draw_transfer_plane(r1, v1, R);
end

% Transfer arc
plt = plot3(R(1,:), R(2,:), R(3,:), 'LineWidth', 2.0);
plt.DisplayName = 'Transfer arc';

% Endpoints
c = get(plt,'Color');
plot3(r1(1), r1(2), r1(3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', c, 'DisplayName','r_1');
plot3(r2(1), r2(2), r2(3), 'x', 'MarkerSize', 10, 'LineWidth', 1.6, 'Color', c, 'DisplayName','r_2');

% Velocity arrows
if opt.ShowVel
    s = quiver_scale(R,V);                 % scale for visibility
    quiver3(r1(1), r1(2), r1(3), v1(1), v1(2), v1(3), s, 'LineWidth',1.4, 'DisplayName','v_1');
    v_end = V(:,end);
    quiver3(R(1,end), R(2,end), R(3,end), v_end(1), v_end(2), v_end(3), s, 'LineWidth',1.4, 'DisplayName','v(t=\Deltat)');
end

legend('Location','bestoutside');
view(35, 20); rotate3d on;
hold off;

end % main function

% --------- helpers ----------
function draw_body(R, nameStr)
    [xs, ys, zs] = sphere(64);
    s = surf(R*xs, R*ys, R*zs, 'EdgeAlpha',0.15, 'FaceAlpha',0.1, ...
             'EdgeColor',[0.4 0.4 0.4], 'FaceColor',[0.5 0.6 0.9], ...
             'DisplayName',char(nameStr));
    uistack(s,'bottom');
    if isempty(nameStr)
        set(s,'DisplayName','Central body');
    end
end

function draw_transfer_plane(r1, v1, R)
    % Plane normal from initial angular momentum
    H = cross(r1, v1);
    n = H / norm(H);
    % Disk radius that comfortably spans the arc
    L = max(range(R,2));
    if ~isfinite(L) || L<=0, L = max(vecnorm(R)); end
    rad = 0.6*max(L, 1);
    % Build orthonormal basis in plane
    % pick any vector not parallel to n
    t = [1;0;0]; if abs(dot(t,n))>0.9, t=[0;1;0]; end
    u = t - n*dot(t,n); u = u/norm(u);
    v = cross(n,u);
    % Parametric disk
    th = linspace(0,2*pi,90);
    rr = linspace(0,rad,12);
    [TH,RR] = meshgrid(th,rr);
    P = u* (RR(:)'.*cos(TH(:)')) + v*(RR(:)'.*sin(TH(:)'));
    X = reshape(P(1,:), size(RR));
    Y = reshape(P(2,:), size(RR));
    Z = reshape(P(3,:), size(RR));
    s = surf(X,Y,Z, 'FaceAlpha',0.05, 'EdgeAlpha',0.05, 'FaceColor',[0 0 0], 'EdgeColor',[0 0 0], 'DisplayName','Transfer plane');
    uistack(s,'bottom');
end

function s = quiver_scale(R,V)
    % pick a scale that keeps arrows readable across scene sizes
    L = max(range(R,2));
    if ~isfinite(L) || L<=0, L = max(vecnorm(R)); end
    Vm = median(vecnorm(V)); if ~isfinite(Vm) || Vm<=0, Vm = 1; end
    % Quiver3 uses the 'scale' factor; we want arrows ≈ few % of scene
    s = 0.05 * L / Vm;
end
