function dx = dymModelFlightSLIPPER(t, x, parms)
% Get parameters
g = parms.g;
% Evalute dimensionless EOMs
yCOMd = x(2);
phid = x(4);

yCOMdd = -g;
phidd = 0;

dx = [yCOMd; yCOMdd;phid ;phidd];
end