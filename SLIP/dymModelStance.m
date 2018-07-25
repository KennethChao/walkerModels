function dx = dymModelStance(t,x,parms)
%
g = parms.g;
k = parms.k;
m = parms.m;
%
l = x(1);
ld = x(2);
theta = x(3);
thetad = x(4);

ldd = l*thetad^2-k/m*(l-1)-g*sin(theta);
thetadd = (-g*l*cos(theta)-2*l*ld*thetad)/l^2;

dx = [ld;ldd;thetad;thetadd];

end