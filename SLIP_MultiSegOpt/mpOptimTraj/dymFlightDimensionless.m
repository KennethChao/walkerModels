function ddx = dymFlightDimensionless(x,dx, parms)
% Get parameters
g = parms.g;
% k = parms.k;
% Get state variables
% x = x(1,:);
% z = x(2,:);
% xd = dx(1,:);
% zd = dx(2,:);
% Evalute dimensionless EOMs of the SLIP model
xdd = 0;
zdd = -g;

ddx = [xdd; zdd];

end