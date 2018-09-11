clc;

addpath('../')

syms g k beta real
 

%% State in Cartesian space at the start of the stance phase
syms l0 l0d l0dd theta0 theta0d theta0dd real

xStance0Polar = [l0;theta0;l0d;theta0d];

[x0Start, z0Start, xd0Start,  zd0Start] = polar2CartesianSLIP(xStance0Polar(1), xStance0Polar(2), xStance0Polar(3), xStance0Polar(4));

xStanceStart = [x0Start; z0Start; xd0Start; zd0Start];

velVecS = [xd0Start; -zd0Start];
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
% syms xF0 zF0 xF0d zF0d xF0dd zF0dd real
% xFlightStart = [xF0; zF0; xF0d; zF0d];

% For legit time

r = roots([-1/2*g, zSEndd, sin(theta0)-zSEnd]);

tFlightEnd = r(2);


variableVector = [lEnd, thetaEnd, lEndd, thetaEndd, g,k, beta];
matlabFunction(tFlightEnd,'File','boundaryConstTFlight','Vars',variableVector);


zFEndd = zSEndd - tFlightEnd*g;
xFEndd = xSEndd;


velVecF = [xFEndd; -zFEndd];
deltaNew = atan2(velVecF(2), velVecF(1));
%% State in Cartesian space at the end of the flight phase
% xStanceEnd = [x(1, parms.phase(1).knotNumber), ...
%     dx(1, parms.phase(1).KnotNumber), ...
%     x(2, parms.phase(1).KnotNumber), ...
%     dx(2, parms.phase(1).KnotNumber)];
% syms xFEnd zFEnd xFEndd zFEndd xFEnddd zFEnddd real
% xFlightEnd = [xFEnd;zFEnd;xFEndd;zFEndd];


% [xF0, xdF0, zF0, zdF0] = polar2CartesianSLIP(xStanceEnd(1), xStanceEnd(2), xStanceEnd(3), xStanceEnd(4));

%% State in Cartesian space at the end of the flight phase

%boundary Condition
cBoundary = simplify([ deltaOld,...
                       ((velVecS(1)^2+velVecS(2)^2) - 1)^2 ]); %v0 = 1                    
                    
                    
c = cBoundary;

variableVector = [l0, theta0, l0d, theta0d, lEnd, thetaEnd, lEndd, thetaEndd, g,k, beta];
matlabFunction(c,'File','boundaryConst','Vars',variableVector);

%% the start of stance phase

x = [l0;theta0];
dx = [l0d;theta0d];
ddx = [l0dd;theta0dd];

xStack = [x;dx;ddx];

jacobianStance0 = simplify(jacobian(c,xStack));
jacobianStance0Pattern = jacobianStance0;
jacobianStance0Pattern(jacobianStance0Pattern~=0) = 1;


% variableVector = [l0, theta0, l0d, theta0d, g,k, beta];
matlabFunction(jacobianStance0,'File','jacobianBoundaryConstX0','Vars',variableVector);
matlabFunction(jacobianStance0Pattern,'File','jacobianBoundaryConstX0Pattern');


