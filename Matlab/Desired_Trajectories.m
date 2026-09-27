%paths
%% circle
function [xd, yd] = Desired_Trajectory(u)
r = 0.5; % Fixed radius (within reach)
theta = (2*pi/6)*u; % fixed speed rotation
xc = 0; yc = 0; % offset from base
xd = xc + r*cos(theta);
yd = yc + r*sin(theta);
end
%% Ellipse
function [xd, yd] = Desired_Trajectory(u)
a = 0.45; % major axis (X direction)
b = 0.25; % minor axis (Y direction)
theta = (2*pi/6)*u; % rotation parameter
xc = 0; yc = 0; % origin-centered
xd = xc + a*cos(theta);
yd = yc + b*sin(theta);
end
%% figure 8
function [xd, yd] = Desired_Trajectory(u)
a = 0.25; % scaling (safe for workspace)
theta = (2*pi/6)*u; % loop frequency
xc = 0; yc = 0;
xd = xc + a*sin(theta);
yd = yc + a*sin(theta).*cos(theta);
end
%% square
function [xd, yd] = Desired_Trajectory(u)
r = 0.35; % half side length
xc = 0; yc = 0;
theta = mod(u,4);

if theta < 1
    xd = xc + r;
    yd = yc + r * (theta);
elseif theta < 2
    xd = xc + r - r * (theta-1);
    yd = yc + r;
elseif theta < 3
    xd = xc - r;
    yd = yc + r - r * (theta-2);
else
    xd = xc - r + r * (theta-3);
    yd = yc - r;
end
end
