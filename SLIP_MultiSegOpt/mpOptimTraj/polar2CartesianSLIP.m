function [x, z, xd, zd] = polar2CartesianSLIP(l, theta, ld, thetad)
%POLAR2CARTESIANSLIP Summary of this function goes here
%   Detailed explanation goes here
z = l .* sin(theta);
zd = ld .* sin(theta) + l .* thetad .* cos(theta);
x = -l * cos(theta);
xd = -ld * cos(theta) + l .* thetad .* sin(theta);
end
