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

% For legit time

r = roots([-1/2*g, zSEndd, sin(theta0)-zSEnd]);

tFlightEnd = r(2);


variableVector = [lEnd, thetaEnd, lEndd, thetaEndd, theta0, g];
matlabFunction(tFlightEnd,'File','getFlightTime','Vars',variableVector);


zFEndd = zSEndd - tFlightEnd*g;
xFEndd = xSEndd;


velVecF = [xFEndd; zFEndd];
deltaNew = atan2(velVecF(2), velVecF(1));

%% State in Cartesian space at the end of the flight phase

%boundary Condition
cBoundary = simplify([  ...%(lEnd -1)^2,...
                        ...((velVecS(1)^2+velVecS(2)^2) - 1)^2,
                        (velVecF'-velVecS').^2
                        ...((velVecS(1)^2+velVecS(2)^2) -(velVecF(1)^2+velVecF(2)^2))^2,...
                        ...deltaNew-deltaOld]); %v0 = 1                    
                        ]);
                    
                    
c = cBoundary

variableVector = [l0, theta0, l0d, theta0d, lEnd, thetaEnd, lEndd, thetaEndd, g,k, beta];
matlabFunction(c,'File','boundaryCost','Vars',variableVector);

%% the start of stance phase

x = [l0;theta0];
dx = [l0d;theta0d];
ddx = [l0dd;theta0dd];

xStack = [x;dx;ddx];

jacobianStance0 = simplify(jacobian(c,xStack));
jacobianStance0Pattern = jacobianStance0;
jacobianStance0Pattern(jacobianStance0Pattern~=0) = 1;


% variableVector = [l0, theta0, l0d, theta0d, g,k, beta];
matlabFunction(jacobianStance0,'File','jacobianBoundaryCostX0','Vars',variableVector);
matlabFunction(jacobianStance0Pattern,'File','jacobianBoundaryCostX0Pattern');


%% the end of stance phase
x = [lEnd;thetaEnd];
dx = [lEndd;thetaEndd];
ddx = [lEnddd;thetaEnddd];

xStack = [x;dx;ddx];

jacobianStanceEnd = simplify(jacobian(c,xStack));
jacobianStanceEndPattern = jacobianStanceEnd;
jacobianStanceEndPattern(jacobianStanceEndPattern~=0) = 1;

% variableVector = [lEnd, thetaEnd, lEndd, thetaEndd,  g,k, beta];
matlabFunction(jacobianStanceEnd,'File','jacobianBoundaryCostXEnd','Vars',variableVector);
matlabFunction(jacobianStanceEndPattern,'File','jacobianBoundaryCostXEndPattern');



%%
velVecFDummy = 1*[xSEndd; zSEndd];

%boundary Condition
cBoundaryFake = simplify([  ...%(lEnd -1)^2,...
                            ...((velVecS(1)^2+velVecS(2)^2) - 1)^2,...
                            (velVecFDummy'-velVecS').^2]); %v0 = 1                    
                    
                    
c = cBoundaryFake

% variableVector = [l0, theta0, l0d, theta0d, lEnd, thetaEnd, g,k, beta];
matlabFunction(c,'File','boundaryCostDummy','Vars',variableVector);

%% the start of stance phase

x = [l0;theta0];
dx = [l0d;theta0d];
ddx = [l0dd;theta0dd];

xStack = [x;dx;ddx];

jacobianStance0 = simplify(jacobian(c,xStack));
% jacobianStance0Pattern = jacobianStance0;
% jacobianStance0Pattern(jacobianStance0Pattern~=0) = 1;


% variableVector = [l0, theta0, l0d, theta0d,   g,k, beta];
matlabFunction(jacobianStance0,'File','jacobianBoundaryCostX0Dummy','Vars',variableVector);
% matlabFunction(jacobianStance0Pattern,'File','jacobianBoundaryConstXS0SinglePatternFake');

%% the end of stance phase
x = [lEnd;thetaEnd];
dx = [lEndd;thetaEndd];
ddx = [lEnddd;thetaEnddd];

xStack = [x;dx;ddx];

jacobianStanceEnd = simplify(jacobian(c,xStack));
% jacobianStanceEndPattern = jacobianStanceEnd;
% jacobianStanceEndPattern(jacobianStanceEndPattern~=0) = 1;

% variableVector = [lEnd, thetaEnd, lEndd, thetaEndd,  g,k, beta];
matlabFunction(jacobianStanceEnd,'File','jacobianBoundaryCostXEndDummy','Vars',variableVector);
% matlabFunction(jacobianStanceEndPattern,'File','jacobianBoundaryCostXEndDummy');

