clc;

addpath('../')

syms h g k real
 

%% State in Cartesian space at the start of the stance phase
syms l ld ldd theta thetad thetadd real

xStancePolar = [l;theta;ld;thetad];

[xStart, zStart, xdStart,  zdStart] = polar2CartesianSLIP(xStancePolar(1), xStancePolar(2), xStancePolar(3), xStancePolar(4));

costStance = simplify((1/2*(xdStart^2 + zdStart^2) - g*zStart - 1/2*k*(l-1)^2 )*h); % 


variableVector = [l, theta, ld, thetad, h, g, k];
matlabFunction(costStance,'File','costStance','Vars',variableVector);


x = [l;theta];
dx = [ld;thetad];
ddx = [ldd;thetadd];

xStack = [x;dx;ddx];

jacobianCostStanceX = simplify(jacobian(costStance,xStack));
variableVector = [l, theta, ld, thetad, h, g, k];
matlabFunction(jacobianCostStanceX,'File','jacobianCostStanceX','Vars',variableVector);

jacobianCostStanceH = simplify(jacobian(costStance,h));
variableVector = [l, theta, ld, thetad, h, g, k];
matlabFunction(jacobianCostStanceH,'File','jacobianCostStanceH','Vars',variableVector);
