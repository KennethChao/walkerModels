clc;

addpath('../')

syms h real
 

%% State in Cartesian space at the start of the stance phase
syms l ld ldd theta thetad thetadd real

xStancePolar = [l;theta;ld;thetad];

[~, ~, xdStart,  zdStart] = polar2CartesianSLIP(xStancePolar(1), xStancePolar(2), xStancePolar(3), xStancePolar(4));

costStance = simplify((xdStart^2 + zdStart^2)*h);


variableVector = [l, theta, ld, thetad, h];
matlabFunction(costStance,'File','costStance','Vars',variableVector);


x = [l;theta];
dx = [ld;thetad];
ddx = [ldd;thetadd];

xStack = [x;dx;ddx];

jacobianCostStanceX = simplify(jacobian(costStance,xStack));
variableVector = [l, theta, ld, thetad, h];
matlabFunction(jacobianCostStanceX,'File','jacobianCostStanceX','Vars',variableVector);

jacobianCostStanceH = simplify(jacobian(costStance,h));
variableVector = [l, theta, ld, thetad, h];
matlabFunction(jacobianCostStanceH,'File','jacobianCostStanceH','Vars',variableVector);

%% State in Cartesian space at the start of the flight phase
% xStanceEnd = [x(1, parms.phase(1).knotNumber), ...
%     dx(1, parms.phase(1).KnotNumber), ...
%     x(2, parms.phase(1).KnotNumber), ...
%     dx(2, parms.phase(1).KnotNumber)];
syms xF zF xFd zFd xFdd zFdd real
xFlightStart = [xF; zF; xFd; zFd];
costFlight = simplify((xFd^2 + zFd^2)*h);

variableVector = [xF, zF, xFd, zFd, h];
matlabFunction(costFlight,'File','costFlight','Vars',variableVector);

x = [xF;zF];
dx = [xFd;zFd];
ddx = [xFdd;zFdd];

xStack = [x;dx;ddx];

jacobianCostFlightX = simplify(jacobian(costFlight,xStack));
variableVector = [xF, zF, xFd, zFd, h];
matlabFunction(jacobianCostFlightX,'File','jacobianCostFlightX','Vars',variableVector);

jacobianCostFlightH = simplify(jacobian(costFlight,h));
variableVector = [xF, zF, xFd, zFd, h];
matlabFunction(jacobianCostFlightH,'File','jacobianCostFlightH','Vars',variableVector);
% 
% %% State in Cartesian space at the end of the flight phase
% % xStanceEnd = [x(1, parms.phase(1).knotNumber), ...
% %     dx(1, parms.phase(1).KnotNumber), ...
% %     x(2, parms.phase(1).KnotNumber), ...
% %     dx(2, parms.phase(1).KnotNumber)];
% syms xFEnd zFEnd xFEndd zFEndd xFEnddd zFEnddd real
% xFlightEnd = [xFEnd;zFEnd;xFEndd;zFEndd];
% 
% 
% % [xF0, xdF0, zF0, zdF0] = polar2CartesianSLIP(xStanceEnd(1), xStanceEnd(2), xStanceEnd(3), xStanceEnd(4));
% 
% %% State in Cartesian space at the end of the flight phase
% % xFlightEnd = [x(1, parms.totalKnotNumber), ...
% %     dx(1, parms.totalKnotNumber), ...
% %     x(2, parms.totalKnotNumber), ...
% %     dx(2, parms.totalKnotNumber)];
% % 
% cTransition = [xFlightEnd - xStanceStart;...
%                  xFlightStart - xStanceEnd];
%              
% %boundary Condition
% cBoundary = simplify((velVecS(1)^2+velVecS(2)^2) - 1); ... %v0 = 1
%     
% c = [cTransition; cBoundary]'
% 
% variableVector = [l0, theta0, l0d, theta0d, lEnd, thetaEnd, lEndd, thetaEndd, xF0, zF0, xF0d, zF0d, xFEnd, zFEnd, xFEndd, zFEndd, g,k, beta];
% matlabFunction(c,'File','boundaryConst','Vars',variableVector);

% %% the start of stance phase
% 
% x = [l0;theta0];
% dx = [l0d;theta0d];
% ddx = [l0dd;theta0dd];
% 
% xStack = [x;dx;ddx];
% 
% jacobianCostStance = simplify(jacobian(c,xStack));
% jacobianStance0Pattern = jacobianCostStance;
% jacobianStance0Pattern(jacobianStance0Pattern~=0) = 1;
% 
% 
% variableVector = [l0, theta0, l0d, theta0d, g,k, beta];
% matlabFunction(jacobianCostStance,'File','jacobianBoundaryConstXS0','Vars',variableVector);
% matlabFunction(jacobianStance0Pattern,'File','jacobianBoundaryConstXS0Pattern');
% 
% %% the end of stance phase
% x = [lEnd;thetaEnd];
% dx = [lEndd;thetaEndd];
% ddx = [lEnddd;thetaEnddd];
% 
% xStack = [x;dx;ddx];
% 
% jacobianStanceEnd = simplify(jacobian(c,xStack));
% jacobianStanceEndPattern = jacobianStanceEnd;
% jacobianStanceEndPattern(jacobianStanceEndPattern~=0) = 1;
% 
% variableVector = [lEnd, thetaEnd, lEndd, thetaEndd, g,k, beta];
% matlabFunction(jacobianStanceEnd,'File','jacobianBoundaryConstXSEnd','Vars',variableVector);
% matlabFunction(jacobianStanceEndPattern,'File','jacobianBoundaryConstXSEndPattern');
% 
% %% the start of flight phase
% x = [xF0;zF0];
% dx = [xF0d;zF0d];
% ddx = [xF0dd;zF0dd];
% 
% xStack = [x;dx;ddx];
% 
% jacobianFlight0 = simplify(jacobian(c,xStack));
% jacobianFlight0Pattern = jacobianFlight0;
% jacobianFlight0Pattern(jacobianFlight0Pattern~=0) = 1;
% 
% variableVector = [xF0, zF0, xF0d, zF0d, g,k, beta];
% matlabFunction(jacobianFlight0,'File','jacobianBoundaryConstXF0','Vars',variableVector);
% matlabFunction(jacobianFlight0Pattern,'File','jacobianBoundaryConstXF0Pattern');
% 
% %% the end of flight phase
% x = [xFEnd;zFEnd];
% dx = [xFEndd;zFEndd];
% ddx = [xFEnddd;zFEnddd];
% 
% xStack = [x;dx;ddx];
% 
% jacobianFlightEnd = simplify(jacobian(c,xStack));
% jacobianFlightEndPattern = jacobianFlightEnd;
% jacobianFlightEndPattern(jacobianFlightEndPattern~=0) = 1;
% 
% variableVector = [xFEnd, zFEnd, xFEndd, zFEndd, g,k, beta];
% matlabFunction(jacobianFlightEnd,'File','jacobianBoundaryConstXFEnd','Vars',variableVector);
% matlabFunction(jacobianFlightEndPattern,'File','jacobianBoundaryConstXFEndPattern');
