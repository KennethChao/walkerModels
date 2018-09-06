addpath('../')

syms l ld ldd theta thetad thetadd u g k

x = [l;theta];
dx = [ld;thetad];
ddx = [ldd;thetadd];


xStack = [x;dx;ddx;u]

ddxOde = dymControlledStanceDimensionless(x,dx, u, g, k)

jacobianStance = simplify(jacobian((ddx-ddxOde),xStack));   


variableVector = [l, theta, ld, thetad, u, g,k];
matlabFunction(jacobianStance,'File','jacobianControlledStanceDym','Vars',variableVector);




% syms xb xbd xbdd zb zbd zbdd g k
% 
% x = [xb;zb];
% dx = [xbd;zbd];
% ddx = [xbdd;zbdd];
% 
% 
% xStack = [x;dx;ddx]
% 
% ddxOde = dymFlightDimensionless(x,dx, g, k)
% 
% jacobianFlight = simplify(jacobian((ddx-ddxOde),xStack));
% 
% variableVector = [xb zb, xbd, zbd, g, k];
% matlabFunction(jacobianFlight,'File','jacobianFlightDym','Vars',variableVector);