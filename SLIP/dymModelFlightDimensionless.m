function dx = dymModelFlightDimensionless(t, x, parms)
% Get parameters
g = parms.g;
% Evalute dimensionless EOMs
yd = x(2);
ydd = -g;

dx = [yd; ydd];

end