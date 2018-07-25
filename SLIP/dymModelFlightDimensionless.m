
function dx = dymModelFlight(t,x,parms)
g = parms.g;
yd = x(2);
ydd = -g;

dx = [yd;ydd];

end