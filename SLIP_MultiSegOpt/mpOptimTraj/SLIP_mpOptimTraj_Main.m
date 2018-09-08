%MAIN.m  --  simple walker trajectory optimization
%
% This script sets up a trajectory optimization problem for a simple model
% of walking, and solves it using OptimTraj. The walking model is a double
% pendulum, with point feet, no ankle torques, impulsive heel-strike (but
% not push-off), and continuous hip torque. Both legs have inertia. Cost
% function is minimize integral of torque-squared.
%
%

clc;
clear;
close all;
% addpath ../../
addpath('./const')
addpath('./helperFunctions')
% ToDO:


% Cartesin2BetaDelta: Skipped
% Gradient:
% - cnst (Periodic)(low priority)
% - cost (low priority)
% Improve Initial Guess


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                  Parameters for the dynamics function                   %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%% Parameter Set
parms.g = 0.66;
parms.beta = 72 / 180 * pi;
parms.k = 20;

parms.delta = 0.1392;


parms.weightBoundary = 1e1;
%% Sys
parms.ndof = 2;
parms.nVarSeg = parms.ndof * 3;
parms.nBoundaryConst = 1;

%% Opt
parms.phase(1).knotNumber = 5;
% parms.phase(2).knotNumber = 21;

totalKnotNumber = 0;
totaHSMCnstNumber = 0;

for i = 1:length(parms.phase)
    totalKnotNumber = totalKnotNumber + parms.phase(i).knotNumber;
    totaHSMCnstNumber = totaHSMCnstNumber + (parms.phase(i).knotNumber-1)/2;
end
    
index = 1;

for i = 1:length(parms.phase)
    parms.phase(i).x0knotNumber = index;
    index = index + parms.phase(i).knotNumber;
end

parms.totalKnotNumber = totalKnotNumber;
parms.totaHSMCnstNumber = totaHSMCnstNumber;
parms.phaseNum = length(parms.phase);
parms.totalVarNumber = parms.totalKnotNumber * parms.nVarSeg + parms.phaseNum;

% Dym
parms.phase(1).dymFunc = @dymStanceDimensionless;

parms.phase(1).jacobianDymFunc = @jacobianStanceDym;

% Cost
parms.phase(1).costFunc = @costStance;
parms.phase(1).jacobianCostX = @jacobianCostStanceX;
parms.phase(1).jacobianCostH = @jacobianCostStanceH;

% Boundary Constraints

parms.phase(1).jacobianBoundaryCostX0 = @jacobianBoundaryCostX0;
parms.phase(1).jacobianBoundaryCostXEnd = @jacobianBoundaryCostXEnd;

parms.phase(1).jacobianBoundaryCostX0Dummy = @jacobianBoundaryCostX0Dummy;
parms.phase(1).jacobianBoundaryCostXEndDummy = @jacobianBoundaryCostXEndDummy;


parms.phase(1).jacobianBoundaryConstX0 = @jacobianBoundaryConstX0;
parms.phase(1).jacobianBoundaryConstX0Pattern = @jacobianBoundaryConstX0Pattern;



% Bounds
parms.phase(1).xlb = [0, parms.beta];
parms.phase(1).dxlb = [-inf, -inf];
parms.phase(1).ddxlb = [-inf, -inf];
parms.phase(1).hlb = 1e-2;

parms.phase(1).xub = [1, pi];
parms.phase(1).dxub = [inf, inf];
parms.phase(1).ddxub = [inf, inf];
parms.phase(1).hub = 10;


% Cechking input
for i = 1:length(parms.phase)
    if mod(parms.phase(i).knotNumber,2)==0
        error('knotNumber should be an odd number')
    end
end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                       Set up function handles                           %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
cost = @(x)costFunction(x, parms); %chk
gCost = @(x)gcostFunction(x, parms); %chk
%
cnst = @(x)constAll(x, parms);
gCnst = @(x)gconstAll(x, parms);
%
% xVec = initialGuess(parms);
cnstPattern = @()GP(parms);

% Assign Function Handle for IPOPT
funcs = {};
funcs.objective = cost;
funcs.gradient = gCost;
funcs.constraints = cnst;
funcs.jacobian = gCnst;
funcs.jacobianstructure =cnstPattern;


% problem.func.dynamics = @(t,x,u)( dymModelStanceDimensionless(t, x,u,parms) );
%
% problem.func.pathObj = @(t,x,u)( cost(t,x) );
%
% problem.func.bndCst = @(t0,x0,tF,xF)( periodicGait(xF,x0,parms) );

%%

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%               Set up bounds on time, state, and control                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
[lb, ub] = inputBounds(parms);
[clb, cub] = constBounds(parms);

options.lb = lb;
options.ub = ub;
options.cl = clb;
options.cu = cub;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%              Create an initial guess for the trajectory                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
xVec = initialGuess(parms);
[x, dx, ddx, h] = extractState(xVec, parms);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           Options:                                      %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
options.ipopt.hessian_approximation = 'limited-memory';
options.ipopt.tol = 1e-4;
options.ipopt.mu_strategy = 'adaptive';

options.ipopt.print_info_string = 'yes';
% options.ipopt.linear_solver = 'ma57';
% options.ipopt.honor_original_bounds = 'no';
% options.ipopt.derivative_test       = 'first-order';
options.ipopt.max_iter = 1500;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           Solve!                                        %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%         if optPars.flag_IPOPT==1
xVec = initialGuess(parms);
[x_Flat2, ~] = ipopt(xVec,funcs,options);

[x, dx, ddx, h] = extractState(x_Flat2, parms);

figure()
plot(x','DisplayName','x')
figure()
plot(dx','DisplayName','dx')
figure()
plot(ddx','DisplayName','ddx')


xStance0Polar = [x(:,1);dx(:,1)];

[x0Start, z0Start, xd0Start,  zd0Start] = polar2CartesianSLIP(xStance0Polar(1), xStance0Polar(2), xStance0Polar(3), xStance0Polar(4));
    
    velVec = [xd0Start, -zd0Start];
    deltaNew = atan2(velVec(2), velVec(1))
    

%         else
%%%%% THE KEY LINE:
% soln = optimTraj(problem);
%
% % Transcription Grid points:
% t = soln(end).grid.time;
% q1 = soln(end).grid.state(1,:);
% q2 = soln(end).grid.state(3,:);
% dq1 = soln(end).grid.state(2,:);
% dq2 = soln(end).grid.state(4,:);
% u = soln(end).grid.control;
%
% % Interpolated solution:
% tInt = linspace(t(1),t(end),10*length(t)+1);
% xInt = soln(end).interp.state(tInt);
% q1Int = xInt(1,:);
% q2Int = xInt(3,:);
% dq1Int = xInt(2,:);
% dq2Int = xInt(4,:);
% uInt = soln(end).interp.control(tInt);
%
% xS = [q1Int(1),dq1Int(1),q2Int(1),q2Int(1)];
%     yS0 = xS( 1) * sin(xS( 3));
%     ySd0 = xS( 2) * sin(xS( 3)) + xS( 1) * xS( 4) * cos(xS( 3));
%     xS0 = -xS( 1) * cos(xS( 3)); %#ok<NASGU>
%     xSd0 = -xS( 2) * cos(xS( 3)) + xS( 1) * xS( 4) * sin(xS( 3));
%
%     velVec = [xSd0, -ySd0];
%     deltaOld = atan2(velVec(2), velVec(1));
%
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% %                     Plot the solution                                   %
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%
% figure(100); clf;
%
% subplot(3,1,1); hold on;
% plot(tInt,q1Int,'r-'); plot(tInt,q2Int,'b-');
% plot([t(1),t(end)],[0,0],'k--','LineWidth',1);
% plot(t,q1,'ro'); plot(t,q2,'bo');
% legend('leg one','leg two')
% xlabel('time (sec)')
% ylabel('angle (rad)')
% title('Leg Angles')
%
% subplot(3,1,2); hold on;
% plot(tInt,dq1Int,'r-'); plot(tInt,dq2Int,'b-');
% plot(t,dq1,'ro'); plot(t,dq2,'bo');
% legend('leg one','leg two')
% xlabel('time (sec)')
% ylabel('rate (rad/sec)')
% title('Leg Angle Rates')
%
% subplot(3,1,3); hold on;
% plot(t,u,'mo'); plot(tInt,uInt,'m-');
% xlabel('time (sec)')
% ylabel('torque (Nm)')
% title('Hip Torque')
function c = costFunction(xVec, parms)
[x, dx, ddx, h] = extractState(xVec, parms);
%
c = costFun(x,dx,h,parms);
end %function end

function g = gcostFunction(xVec, parms)
[x, dx, ddx, h] = extractState(xVec, parms);
%
g = gcostFun(x,dx,h,parms);
end %function end




function [lb, ub] = inputBounds(parms)
lb = ones(1, parms.totalVarNumber) * -inf;
ub = ones(1, parms.totalVarNumber) * inf;

%%
shiftIndex = 0;
for i = 1:parms.phaseNum
    for j = 1:parms.phase(i).knotNumber
        lb(1, (1:parms.ndof)+(j - 1)*parms.nVarSeg+shiftIndex) = parms.phase(i).xlb;
        lb(1, (1:parms.ndof)+(j - 1)*parms.nVarSeg+parms.ndof+shiftIndex) = parms.phase(i).dxlb;
        lb(1, (1:parms.ndof)+(j - 1)*parms.nVarSeg+parms.ndof*2+shiftIndex) = parms.phase(i).ddxlb;
        
        ub(1, (1:parms.ndof)+(j - 1)*parms.nVarSeg+shiftIndex) = parms.phase(i).xub;
        ub(1, (1:parms.ndof)+(j - 1)*parms.nVarSeg+parms.ndof+shiftIndex) = parms.phase(i).dxub;
        ub(1, (1:parms.ndof)+(j - 1)*parms.nVarSeg+parms.ndof*2+shiftIndex) = parms.phase(i).ddxub;
    end
    shiftIndex = shiftIndex + parms.phase(i).knotNumber * parms.nVarSeg;
end
for i = (1:length(parms.phase))
    lb(1, end-i+1) = parms.phase(i).hlb;
    ub(1, end-i+1) = parms.phase(i).hub;
end

lb(1) = 1;
ub(1) = 1;
% 
lb(2) = parms.beta;
ub(2) = parms.beta;

% lb(3) = -1 * cos(parms.beta-parms.delta);
% ub(3) = -1 * cos(parms.beta-parms.delta);
% % 
% lb(4) = 1 * sin(parms.beta-parms.delta);
% ub(4) = 1 * sin(parms.beta-parms.delta);
end %function end

function [clb, cub] = constBounds(parms)

%%
% Bounds of Kinematic Constraints
nHSM = 2;
relativeDegree = 2;
clbKine = zeros(parms.totaHSMCnstNumber*parms.ndof*nHSM*relativeDegree, 1);
cubKine = zeros(parms.totaHSMCnstNumber*parms.ndof*nHSM*relativeDegree, 1);

% Bounds of Dynamic Constraints
clbDym = zeros((parms.totalKnotNumber)*parms.ndof, 1);
cubDym = zeros((parms.totalKnotNumber)*parms.ndof, 1);

% Bounds of Boundary Constraints
% tol = 1e-2;
% clbBoundary = -tol*ones(parms.nBoundaryConst, 1);
% cubBoundary = tol*ones(parms.nBoundaryConst, 1);
clbBoundary = 0.2;
cubBoundary = inf;
% cubBoundary(end) = inf;

clb = [clbKine; clbDym; clbBoundary];
cub = [cubKine; cubDym; cubBoundary];

% clb = [clbKine; clbDym];
% cub = [cubKine; cubDym];

% clb = clbBoundary;
% cub = cubBoundary;

end %function end

function c = constAll(xVec, parms)
[x, dx, ddx, h] = extractState(xVec, parms);

c0 = constKineHSM(x, dx, ddx, h, parms);
c1 = constDym(x, dx, ddx, parms);
c2 = constBoundary(x, dx, parms);
c = [c0 c1 c2];
% c = [c2;];
end %function end

function g = gconstAll(xVec, parms)
[x, dx, ddx, h] = extractState(xVec, parms);
g0=gconstKineHSM(x, dx, ddx, h, parms);
g1=gconstDym(x,dx,ddx,parms);
g2=gconstBoundary(x,dx,ddx,parms);
g = [g0;g1;g2];%g1;g2;g3
% g = [g2];%g1;g2;g3
end %function end
%
function Pattern = GP(parms)
G0 = gconstKineHSMPattern(parms);
G1 = gconstDymPattern(parms);
G2 = gconstBoundaryPattern(parms);
Pattern = [G0;G1;G2];
% Pattern = G2;
end
function xVec = state2FreeVariableVector(x, dx, ddx, h, parms)

xVec = zeros(parms.totalKnotNumber*parms.nVarSeg+length(parms.phase), 1);
for i = 1:parms.totalKnotNumber
    xSegment = [x(:, i); ...
        dx(:, i); ...
        ddx(:, i); ...
        ];
    xVec((i - 1)*parms.nVarSeg+(1:parms.nVarSeg), 1) = xSegment;
    
end

% xVec(end-1) = h(1);
% xVec(end) = h(1);
for i = 1:length(parms.phase)
     xVec(end-i+1) = h(i);
end

end
function [x, dx, ddx, h] = extractState(aVec, parms)
x = zeros(parms.ndof, parms.totalKnotNumber);
dx = zeros(parms.ndof, parms.totalKnotNumber);
ddx = zeros(parms.ndof, parms.totalKnotNumber);

for i = 1:parms.totalKnotNumber
    x(:, i) = aVec((i - 1)*parms.nVarSeg+(1:parms.ndof), 1);
    dx(:, i) = aVec((i - 1)*parms.nVarSeg+(1:parms.ndof)+parms.ndof, 1);
    ddx(:, i) = aVec((i - 1)*parms.nVarSeg+(1:parms.ndof)+parms.ndof*2, 1);
end

for i = 1:length(parms.phase)
     h(i) = aVec(end-i+1);
end

end

function xVec = initialGuess(parms)
% h
% h = [0.2, 0.2];
h = [0.2];
x1 = [linspace(1, 1, parms.phase(1).knotNumber); ...
    linspace(parms.beta, parms.beta*2, parms.phase(1).knotNumber);];


% t2 = linspace(1,h(2),parms.phase(2).knotNumber);

% x2 = [linspace(1, 1, parms.phase(2).knotNumber); ...
%     linspace(parms.beta, parms.beta*2, parms.phase(2).knotNumber);];
x2 = [];
x = [x1, x2];

% dx
dx = [0.1*ones(1, parms.totalKnotNumber);
        -0.1*ones(1, parms.totalKnotNumber)
    ];
% ddx
ddx = [0.0*ones(1, parms.totalKnotNumber);
        0.00*ones(1, parms.totalKnotNumber)
    ];


xVec = state2FreeVariableVector(x, dx, ddx, h, parms);

end