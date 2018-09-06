function ddx = dymControlledStanceDimensionless(x,dx, g, u, k)
% Get parameters
% g = parms.g;
% k = parms.k;
% Get state variables
l = x(1,:);
theta = x(2,:);
ld = dx(1,:);
thetad = dx(2,:);
% Evalute dimensionless EOMs of the SLIP model
ldd = l .* thetad.^2 - k * (l - 1) - g * sin(theta) + u;
thetadd = (-g * l .* cos(theta) - 2 * l .* ld .* thetad) ./ (l.^2);

ddx = [ldd; thetadd];

end