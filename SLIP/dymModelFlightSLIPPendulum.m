function dx = dymModelFlightSLIPPendulum(t, x, parms)
% Get parameters
g = parms.g;
% Evalute dimensionless EOMs
yd = x(2);
phid = x(4);

ydd = -g;
phidd = 0;

dx = [yd; ydd;phid ;phidd];
end