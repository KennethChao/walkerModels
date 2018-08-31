addpath('../')

syms l ld ldd theta thetad thetadd g k

x = [l;theta];
dx = [ld;thetad];
ddx = [ldd;thetadd];


xStack = [x;dx;ddx]

ddxOde = dymStanceDimensionless(x,dx, g, k)

jacobianStance = simplify(jacobian((ddx-ddxOde),xStack));   


variableVector = [l, theta, ld, thetad, g,k];
matlabFunction(jacobianStance,'File','jacobianStanceDym','Vars',variableVector);

xStack = [x;dx;ddx]


syms xb xbd xbdd zb zbd zbdd g k

x = [xb;zb];
dx = [xbd;zbd];
ddx = [xbdd;zbdd];


xStack = [x;dx;ddx]

ddxOde = dymFlightDimensionless(x,dx, g, k)

jacobianFlight = simplify(jacobian((ddx-ddxOde),xStack));

variableVector = [xb zb, xbd, zbd, g, k];
matlabFunction(jacobianFlight,'File','jacobianFlightDym','Vars',variableVector);