clc;

addpath('../')

syms g k beta real
 

%% State in Cartesian space at the start of the stance phase
syms l0 l0d l0dd theta0 theta0d theta0dd real

xStance0Polar = [l0;theta0;l0d;theta0d];

[x0Start, z0Start, xd0Start,  zd0Start] = polar2CartesianSLIP(xStance0Polar(1), xStance0Polar(2), xStance0Polar(3), xStance0Polar(4));

xStanceStart = [x0Start; z0Start; xd0Start; zd0Start];

velVecS = [xd0Start, -zd0Start];
deltaOld = atan2(velVecS(2), velVecS(1));

%% State in Cartesian space at the end of the stance phase
syms lEnd lEndd lEnddd thetaEnd thetaEndd thetaEnddd real

xStanceEndPolar = [lEnd;thetaEnd;lEndd;thetaEndd];

[xSEnd, zSEnd, xSEndd, zSEndd] = polar2CartesianSLIP(xStanceEndPolar(1), xStanceEndPolar(2), xStanceEndPolar(3), xStanceEndPolar(4));

xStanceEnd = [xSEnd; xSEndd; zSEnd; zSEndd];

%% State in Cartesian space at the start of the flight phase
% xStanceEnd = [x(1, parms.phase(1).knotNumber), ...
%     dx(1, parms.phase(1).KnotNumber), ...
%     x(2, parms.phase(1).KnotNumber), ...
%     dx(2, parms.phase(1).KnotNumber)];
syms xF0 zF0 xF0d zF0d xF0dd zF0dd real
xFlightStart = [xF0; zF0; xF0d; zF0d];


%% State in Cartesian space at the end of the flight phase
% xStanceEnd = [x(1, parms.phase(1).knotNumber), ...
%     dx(1, parms.phase(1).KnotNumber), ...
%     x(2, parms.phase(1).KnotNumber), ...
%     dx(2, parms.phase(1).KnotNumber)];
syms xFEnd zFEnd xFEndd zFEndd xFEnddd zFEnddd real
xFlightEnd = [xFEnd;zFEnd;xFEndd;zFEndd];


% [xF0, xdF0, zF0, zdF0] = polar2CartesianSLIP(xStanceEnd(1), xStanceEnd(2), xStanceEnd(3), xStanceEnd(4));

%% State in Cartesian space at the end of the flight phase
% xFlightEnd = [x(1, parms.totalKnotNumber), ...
%     dx(1, parms.totalKnotNumber), ...
%     x(2, parms.totalKnotNumber), ...
%     dx(2, parms.totalKnotNumber)];
% 
cTransition = [xFlightEnd - xStanceStart;...
                 xFlightStart - xStanceEnd];
             
%boundary Condition
cBoundary = simplify([  xStance0Polar(1)-1;
                        xStance0Polar(2)-beta;
                        xStanceEndPolar(1)-1;
                        (velVecS(1)^2+velVecS(2)^2) - 1;
                        xFEnd-0.0]); %v0 = 1
                        
c = [cTransition; cBoundary]'

variableVector = [l0, theta0, l0d, theta0d, lEnd, thetaEnd, lEndd, thetaEndd, xF0, zF0, xF0d, zF0d, xFEnd, zFEnd, xFEndd, zFEndd, g,k, beta];
matlabFunction(c,'File','boundaryConst','Vars',variableVector);

%% the start of stance phase

x = [l0;theta0];
dx = [l0d;theta0d];
ddx = [l0dd;theta0dd];

xStack = [x;dx;ddx];

jacobianStance0 = simplify(jacobian(c,xStack));
jacobianStance0Pattern = jacobianStance0;
jacobianStance0Pattern(jacobianStance0Pattern~=0) = 1;


variableVector = [l0, theta0, l0d, theta0d, g,k, beta];
matlabFunction(jacobianStance0,'File','jacobianBoundaryConstXS0','Vars',variableVector);
matlabFunction(jacobianStance0Pattern,'File','jacobianBoundaryConstXS0Pattern');

%% the end of stance phase
x = [lEnd;thetaEnd];
dx = [lEndd;thetaEndd];
ddx = [lEnddd;thetaEnddd];

xStack = [x;dx;ddx];

jacobianStanceEnd = simplify(jacobian(c,xStack));
jacobianStanceEndPattern = jacobianStanceEnd;
jacobianStanceEndPattern(jacobianStanceEndPattern~=0) = 1;

variableVector = [lEnd, thetaEnd, lEndd, thetaEndd, g,k, beta];
matlabFunction(jacobianStanceEnd,'File','jacobianBoundaryConstXSEnd','Vars',variableVector);
matlabFunction(jacobianStanceEndPattern,'File','jacobianBoundaryConstXSEndPattern');

%% the start of flight phase
x = [xF0;zF0];
dx = [xF0d;zF0d];
ddx = [xF0dd;zF0dd];

xStack = [x;dx;ddx];

jacobianFlight0 = simplify(jacobian(c,xStack));
jacobianFlight0Pattern = jacobianFlight0;
jacobianFlight0Pattern(jacobianFlight0Pattern~=0) = 1;

variableVector = [xF0, zF0, xF0d, zF0d, g,k, beta];
matlabFunction(jacobianFlight0,'File','jacobianBoundaryConstXF0','Vars',variableVector);
matlabFunction(jacobianFlight0Pattern,'File','jacobianBoundaryConstXF0Pattern');

%% the end of flight phase
x = [xFEnd;zFEnd];
dx = [xFEndd;zFEndd];
ddx = [xFEnddd;zFEnddd];

xStack = [x;dx;ddx];

jacobianFlightEnd = simplify(jacobian(c,xStack));
jacobianFlightEndPattern = jacobianFlightEnd;
jacobianFlightEndPattern(jacobianFlightEndPattern~=0) = 1;

variableVector = [xFEnd, zFEnd, xFEndd, zFEndd, g,k, beta];
matlabFunction(jacobianFlightEnd,'File','jacobianBoundaryConstXFEnd','Vars',variableVector);
matlabFunction(jacobianFlightEndPattern,'File','jacobianBoundaryConstXFEndPattern');
