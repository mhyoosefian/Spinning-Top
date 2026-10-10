clc; clear; close all;
% Simulates a spinning top whose tip rests on the floor, using the modified
% Lagrange method (Euler angles or quaternions), and plays a 3D animation.
% A friction torque about the symmetry axis slowly reduces the spin; when
% the spin is too low, the top tilts further until its rim touches the
% floor. Set saveVideo = 1 to also save the animation as an MP4 file.

global m h Moment g Ixx Iyy Izz
P = topProfile();                         % shape and mass properties
m = P.m; h = P.h; Ixx = P.Ixx; Iyy = P.Ixx; Izz = P.Izz; g = 9.81; Moment = 0;

%_________________________Settings_________________________________________
method = 'euler';        % 'euler' or 'quat'
spin0 = 150;             % initial spin about the symmetry axis (rad/s)
tilt0 = 8;               % initial tilt from the vertical (deg)
friction = 2.4e-6;       % friction torque = -friction * spin (N m s)
Tmax = 30;               % longest simulated time (s)
speed = 1/4;             % playback speed (1/4 = four times slower than real time)
fps = 60;                % frames per second of the animation
saveVideo = 0;           % 1 to save the animation as an MP4 file

%_________________________Simulation_______________________________________
% States: [psi; theta; phi; u] (Euler) or [e0; e1; e2; e3; u] (quaternions),
% where u is the angular velocity in the x2y2z2 (Euler) or body (quaternion)
% frame. The simulation stops when the rim touches the floor.
angles0 = [0; deg2rad(tilt0); 0];
w0 = [0; 0; spin0];
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-10, 'MaxStep', 1e-3, 'Events', @(t, z) floorContact(t, z, method, P));
switch method
    case 'euler'
        sol = ode45(@(t, z) dynEuler(t, z, friction), [0 Tmax], [angles0; w0], options);
    case 'quat'
        % The quaternion functions were generated before the sign of the
        % potential energy was corrected (V = +m*g*h*cos(theta)), so their
        % gravity points up. Flipping g gives the same top as the Euler code.
        g = -g;
        sol = ode45(@(t, z) dynQuat(t, z, friction), [0 Tmax], [rot2quat(rotZYZ(0, angles0(2), 0)); w0], options);
        g = -g;
end
tFall = sol.x(end);
t = (0:speed/fps:tFall).';
z = deval(sol, t).';
N = numel(t); R = zeros(3, 3, N); spin = z(:, end);
for k = 1:N
    if strcmp(method, 'euler')
        R(:, :, k) = rotZYZ(z(k, 1), z(k, 2), z(k, 3));
    else
        R(:, :, k) = quat2rot(z(k, 1:4)/norm(z(k, 1:4)));
    end
end
fell = tFall < Tmax;

% After the rim touches the floor the top rolls on its rim and stops. This
% last part is not simulated: the precession slows down linearly to zero
% while the top rolls without slipping (spin = -precession/cos(theta)).
tEnd = []; REnd = [];
if fell
    az = @(Rk) atan2(Rk(2, 3), Rk(1, 3));
    rate = angdiff(az(R(:, :, end-1)), az(R(:, :, end)))/(t(2) - t(1));
    Te = .4; tau = (speed/fps:speed/fps:Te + .3).';
    dpsi = rate*(min(tau, Te) - min(tau, Te).^2/(2*Te));
    Rz = @(a) [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1];
    REnd = zeros(3, 3, numel(tau));
    for k = 1:numel(tau)
        REnd(:, :, k) = Rz(dpsi(k))*R(:, :, end)*Rz(-dpsi(k)/cos(P.thetaContact));
    end
    tEnd = tFall + tau;
end

%_________________________Animation________________________________________
tAll = [t; tEnd]; RAll = cat(3, R, REnd);
readout = cell(numel(tAll), 1);
for k = 1:numel(tAll)
    th = acosd(RAll(3, 3, k));
    if k <= N
        s = sprintf('spin = %4.1f rev/s', spin(k)/(2*pi));
    else
        s = 'the rim touches the floor';
    end
    readout{k} = sprintf('Spinning top (%g\\times slow motion)\nt = %4.2f s,  \\theta = %4.1f\\circ\n%s', 1/speed, tAll(k), th, s);
end
opts.speed = speed; opts.readout = readout; opts.trail = .4;
if saveVideo, opts.video = sprintf('spinning_top_%s.mp4', method); end
animateTop(tAll, RAll, opts);
if fell
    fprintf('The rim touched the floor at t = %.2f s (spin %.0f rad/s).\n', tFall, spin(end));
else
    fprintf('The top was still spinning at t = %.0f s.\n', Tmax);
end


%_________________________Helper functions_________________________________
function dz = dynEuler(t, z, c)
% ModlagDyn with a friction torque about the symmetry axis.
global Moment
Moment = -c*z(6);
dz = ModlagDyn(t, z);
end

function dz = dynQuat(t, z, c)
% ModlagDyn_quat with a friction torque about the symmetry axis.
global Moment
Moment = -c*z(7);
dz = ModlagDyn_quat(t, z);
end

function [value, isterminal, direction] = floorContact(~, z, method, P)
% Zero when the tilt reaches the angle at which the rim touches the floor.
if strcmp(method, 'euler')
    theta = z(2);
else
    e = z(1:4)/norm(z(1:4));
    theta = acos(max(-1, min(1, e(1)^2 - e(2)^2 - e(3)^2 + e(4)^2)));
end
value = theta - P.thetaContact; isterminal = 1; direction = 1;
end

function R = rotZYZ(psi, theta, phi)
% Rotation from the body frame to the space frame for ZYZ Euler angles.
Rz = @(a) [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1];
Ry = @(a) [cos(a) 0 sin(a); 0 1 0; -sin(a) 0 cos(a)];
R = Rz(psi)*Ry(theta)*Rz(phi);
end

function R = quat2rot(e)
% Body-to-space rotation. The matrix R in Symb_Deriv.m maps space to body,
% so this is its transpose.
e0 = e(1); e1 = e(2); e2 = e(3); e3 = e(4);
Rsb = [e0^2+e1^2-e2^2-e3^2, 2*(e1*e2+e0*e3), 2*(e1*e3-e0*e2);
       2*(e1*e2-e0*e3), e0^2+e2^2-e1^2-e3^2, 2*(e2*e3+e0*e1);
       2*(e1*e3+e0*e2), 2*(e2*e3-e0*e1), e0^2+e3^2-e1^2-e2^2];
R = Rsb.';
end

function e = rot2quat(R)
% Quaternion [e0; e1; e2; e3] of a body-to-space rotation (inverse of quat2rot).
e0 = sqrt(max(0, 1 + trace(R)))/2;
e = [e0; (R(3,2)-R(2,3))/(4*e0); (R(1,3)-R(3,1))/(4*e0); (R(2,1)-R(1,2))/(4*e0)];
end

function d = angdiff(a, b)
% Smallest signed difference b - a between two angles.
d = mod(b - a + pi, 2*pi) - pi;
end
