addpath('../')

syms l ld ldd theta thetad thetadd g k beta
 
% x = [l;theta];
% dx = [ld;thetad];
% ddx = [ldd;thetadd];


%% State in Cartesian space at the start of the stance phase

xStance0Polar = [l;theta;ld;thetad];

[x0Start, xd0Start, zSStart, zd0Start] = polar2CartesianSLIP(xStance0Polar(1), xStance0Polar(2), xStance0Polar(3), xStance0Polar(4));

xStance0 = [x0Start; xd0Start; zSStart; zd0Start];

velVecS = [xd0Start, -zd0Start];
deltaOld = atan2(velVecS(2), velVecS(1));

%% State in Cartesian space at the end of the stance phase
% xStanceEnd = [x(1, parms.phase(1).knotNumber), ...
%     dx(1, parms.phase(1).KnotNumber), ...
%     x(2, parms.phase(1).KnotNumber), ...
%     dx(2, parms.phase(1).KnotNumber)];
syms xEnd zEnd xdEnd zdEnd xddEnd zddEnd
xStanceEnd = [xEnd;zEnd;xdEnd;zdEnd];


% [xF0, xdF0, zF0, zdF0] = polar2CartesianSLIP(xStanceEnd(1), xStanceEnd(2), xStanceEnd(3), xStanceEnd(4));

%% State in Cartesian space at the end of the flight phase
% xFlightEnd = [x(1, parms.totalKnotNumber), ...
%     dx(1, parms.totalKnotNumber), ...
%     x(2, parms.totalKnotNumber), ...
%     dx(2, parms.totalKnotNumber)];
% 
cTransition = xStanceEnd - xStance0;
%boundary Condition
cBoundary = simplify([xStance0Polar(1) - 1; ... %l0=1
    xStance0Polar(2) - beta; ... %theta0=beta
    (velVecS(1)^2+velVecS(2)^2) - 1; ... %v0 = 1
    ]);

c = [cTransition; cBoundary]'
%%

x = [l;theta];
dx = [ld;thetad];
ddx = [ldd;thetadd];

xStack = [x;dx;ddx];

jacobianStance = simplify(jacobian(c,xStack));


variableVector = [l, theta, ld, thetad, g,k, beta];
matlabFunction(jacobianStance,'File','jacobianPeriodicX0','Vars',variableVector);


%%
x = [xEnd;zEnd];
dx = [xdEnd;zdEnd];
ddx = [xddEnd;zddEnd];

xStack = [x;dx;ddx];

jacobianFlight = simplify(jacobian(c,xStack));


variableVector = [xEnd, zEnd, xdEnd, zdEnd, g,k, beta];
matlabFunction(jacobianFlight,'File','jacobianPeriodicXEnd','Vars',variableVector);
