function dx = dymModelStanceDimensionless(t, x, u, parms)
% Get parameters
g = parms.g;
k = parms.k;
% Get state variables
l = x(1,:);
ld = x(2,:);
theta = x(3,:);
thetad = x(4,:);
% Evalute dimensionless EOMs of the SLIP model
ldd = l .* thetad.^2 - k * (l - 1) - g * sin(theta);
thetadd = (-g * l .* cos(theta) - 2 * l .* ld .* thetad) ./ (l.^2);

dx = [ld; ldd; thetad; thetadd];

end