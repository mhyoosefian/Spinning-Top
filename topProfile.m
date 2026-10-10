function P = topProfile()
% TOPPROFILE  Shape of the spinning top and its mass properties.
%   The top is a body of revolution about its symmetry axis (body z axis)
%   with its tip at the origin. Its mass, center of mass and moments of
%   inertia are computed from the shape assuming a uniform density, so the
%   animation and the simulation describe the same object.
%
%   P.z, P.r   profile of the body: radius r (m) at height z (m)
%   P.zs, P.rs stem: heights and radius (m)
%   P.m, P.h   mass (kg) and distance of the center of mass from the tip (m)
%   P.Ixx, P.Izz  moments of inertia about the center of mass (kg m^2)
%   P.thetaContact  tilt angle at which the rim touches the floor (rad)

P.z = [0 .002 .012 .022 .030 .033 .035 .037 .039 .041 .043 .044 .045];
P.r = [0 .002 .010 .020 .028 .0305 .031 .030 .027 .021 .013 .006 .004];
P.zs = [.045 .062]; P.rs = .0035;
rho = 1150;                                   % density of plastic (kg/m^3)

% Body: stack of thin disks
z = linspace(0, P.z(end), 2000); r = interp1(P.z, P.r, z); dz = z(2) - z(1);
dm = rho*pi*r.^2*dz;
mb = sum(dm); zb = sum(dm.*z)/mb;
Izz = sum(dm.*r.^2/2);
Ixx_tip = sum(dm.*(r.^2/4 + z.^2));
% Stem: solid cylinder
ls = diff(P.zs); ms = rho*pi*P.rs^2*ls; zcs = mean(P.zs);
Izz = Izz + ms*P.rs^2/2;
Ixx_tip = Ixx_tip + ms*(P.rs^2/4 + ls^2/12 + zcs^2);

P.m = mb + ms;
P.h = (mb*zb + ms*zcs)/P.m;
P.Izz = Izz;
P.Ixx = Ixx_tip - P.m*P.h^2;                  % about the center of mass
% The rim touches the floor when some point (r, z) of the profile reaches
% the floor: z*cos(theta) = r*sin(theta).
P.thetaContact = min(atan2(P.z(2:end), P.r(2:end)));
end
