function dx = dymModelStanceSLIPPendulum(t, x, parms)
% Get parameters
g = parms.g;
k = parms.k;
mf = parms.mf;
rc = parms.rc;
% Get state variables
l = x(1);
ld = x(2);
theta = x(3);
thetad = x(4);
phi = x(5);
phid = x(6);

M = [1+mf,                0,                   rc*sin(phi-theta);
     0                    (1+mf)*l^2,          -rc*cos(phi-theta);
     rc*sin(phi-theta)    -rc*cos(phi-theta)   rc^2              ];

% Evalute dimensionless EOMs of the SLIP model
b = zeros(3,1);

b(1) = -rc*phid*cos(phi-theta)*(phid-thetad)-(1+mf)*g*sin(theta) + (1+mf)*l*thetad^2-k*(l-1)-rc*phid*thetad*cos(phi-theta);
b(2) = -2*(1+mf)*l*ld*thetad + rc*phid*ld*cos(phi-theta)-rc*phid*l*sin(phi-theta)*(phid-thetad)-(1+mf)*g*l*cos(theta) + rc*phid*(-ld*cos(phi-theta)-l*thetad*sin(phi-theta));
b(3) = -rc*ld*cos(phi-theta)*(phid-thetad)+rc*(ld*thetad*cos(phi-theta)-l*thetad*sin(phi-theta)*(phid-thetad))  ...
       + g*rc*cos(phi)+rc*phid*(ld*cos(phi-theta)+l*thetad*sin(phi-theta));

xdd = M^(-1)*b;

ldd = xdd(1);
thetadd = xdd(2);
phidd = xdd(3);
% phidd = b(3)/rc^2;

dx = [ld; ldd; thetad; thetadd; phid; phidd];
% dx = [0; 0; 0; 0; phid; phidd];

end