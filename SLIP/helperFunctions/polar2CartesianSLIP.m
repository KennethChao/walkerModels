function [x, z, xd, zd] = polar2CartesianSLIP(l, theta, ld, thetad)
%POLAR2CARTESIANSLIP function to convert motion from polar cooridnate to
% Cartesian space
%   Covert the COM motions in terms of (l, ld, theta, thetad) to
%   (x, xd, z, zd).

z = l .* sin(theta);
zd = ld .* sin(theta) + l .* thetad .* cos(theta);
x = -l * cos(theta);
xd = -ld * cos(theta) + l .* thetad .* sin(theta);
end
