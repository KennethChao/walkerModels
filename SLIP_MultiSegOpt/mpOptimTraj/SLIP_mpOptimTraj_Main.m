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
% addpath ../../
addpath('./const')
addpath('./helperFunctions')
% ToDO:
% 08/30 - HSM constraints, pattern, tests
    % constraints refactored, 
    % pattern done
    % test with different knot number done
    % test with fewest knot number done
    



% Check vector flatten or unflatten: Done
% Pola2Cartesin (Vector or Matrix): Done
% Cartesin2BetaDelta: Skipped
% Gradient:
% - cnst (Kine)
% - cnst (Dym)
% - cnst (Periodic)(low priority)
% - cost (low priority)
% Improve Initial Guess


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                  Parameters for the dynamics function                   %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%% Parameter Set
parms.g = 0.46;
parms.beta = 72 / 180 * pi;
parms.k = 16;

parms.delta0 = 0.1;

%% Sys
parms.ndof = 2;
parms.nVarSeg = parms.ndof * 3;
parms.nPeriodicConst = 10;

%% Opt
parms.phase(1).knotNumber = 5;
parms.phase(2).knotNumber = 5;

totalKnotNumber = 0;
totaHSMCnstNumber = 0;

for i = 1:length(parms.phase)
    totalKnotNumber = totalKnotNumber + parms.phase(i).knotNumber;
    totaHSMCnstNumber = totaHSMCnstNumber + (parms.phase(i).knotNumber-1)/2;
end

parms.totalKnotNumber = totalKnotNumber;
parms.totaHSMCnstNumber = totaHSMCnstNumber;
parms.phaseNum = length(parms.phase);
parms.totalVarNumber = parms.totalKnotNumber * parms.nVarSeg + parms.phaseNum;

% Dym
parms.phase(1).dymFunc = @dymStanceDimensionless;
parms.phase(2).dymFunc = @dymFlightDimensionless;

parms.phase(1).jacobianDymFunc = @jacobianStanceDym;
parms.phase(2).jacobianDymFunc = @jacobianFlightDym;

% Bounds
parms.phase(1).xlb = [0, parms.beta];
parms.phase(1).dxlb = [-inf, -inf];
parms.phase(1).ddxlb = [-inf, -inf];
parms.phase(1).hlb = 1e-2;

parms.phase(1).xub = [inf, pi];
parms.phase(1).dxub = [inf, inf];
parms.phase(1).ddxub = [inf, inf];
parms.phase(1).hub = 10;

parms.phase(2).xlb = [-inf, 0];
parms.phase(2).dxlb = [-inf, -inf];
parms.phase(2).ddxlb = [-inf, -inf];
parms.phase(2).hlb = 1e-2;

parms.phase(2).xub = [inf, inf];
parms.phase(2).dxub = [inf, inf];
parms.phase(2).ddxub = [inf, inf];
parms.phase(2).hub = 10;

% Cechking input
for i = 1:length(parms.phase)
    if mod(parms.phase(i).knotNumber,2)==0
        error('knotNumber should be an odd number')
    end
end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                       Set up function handles                           %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
cost = @(x)costTest(x); %chk
gCost = @(x)gcostTest(x); %chk
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
options.ipopt.derivative_test       = 'first-order';
options.ipopt.max_iter = 1;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           Solve!                                        %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%         if optPars.flag_IPOPT==1
xVec = initialGuess(parms);
[x_Flat2, ~] = ipopt(xVec,funcs,options);
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
c = costFun(h, dx);
end %function end

function c = constAll(xVec, parms)
[x, dx, ddx, h] = extractState(xVec, parms);

% c0 = constKineHSM(x, dx, ddx, h, parms);
c1 = constDym(x, dx, ddx, parms);
% c2 = constPeriodic(x, dx, parms);
% c = [c0; c1; c2];
c = [c1;];
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


end %function end

function [clb, cub] = constBounds(parms)

%%
% Bounds of Kinematic Constraints
nHSM = 2;
relativeDegree = 2;
clbKine = zeros(parms.totaHSMCnstNumber*parms.ndof*nHSM*relativeDegree, 1);
cubKine = zeros(parms.totaHSMCnstNumber*parms.ndof*nHSM*relativeDegree, 1);

% Bounds of Dynamic Constraints
clbDym = zeros(parms.totalKnotNumber*parms.ndof, 1);
cubDym = zeros(parms.totalKnotNumber*parms.ndof, 1);

% Bounds of Periodic Constraints
clbPeriodic = zeros(parms.nPeriodicConst, 1);
cubPeriodic = zeros(parms.nPeriodicConst, 1);

% clb = [clbKine; clbDym; clbPeriodic];
% cub = [cubKine; cubDym; cubPeriodic];

clb = [clbDym]';
cub = [cubDym]';
end %function end


function g = gconstAll(xVec, parms)
[x, dx, ddx, h] = extractState(xVec, parms);
% g0=gconstKineHSM(x, dx, ddx, h, parms);
g1=gconstDym(x,dx,ddx,parms);
% g2=gContact(aVec, aInd, cfg);
% g3=gPeriodic(aVec, aInd, cfg);
% g = [g0;g1;g2;g3];%g1;g2;g3
g = [g1];%g1;g2;g3
end %function end
%
function g0 = GP(parms)
% xVec = initialGuess(parms);
% [x, dx, ddx, h] = extractState(xVec, parms);
% g0=gconstKineHSM(x, dx, ddx, h, parms, true);
% g0 = gconstPatternKineHSM(parms);
g0 = sparse(ones(20,62));
% g0 = gconstDymPattern(parms);
% g0 = sparse(ones(24,50));
% G0 = gKineHSMPat(cfg);
% G1 = gDymHSMPat(cfg);
% G2 = gContactPat(cfg);
% G3 = gPeriodicPat(cfg);
% Pattern = [G0;G1;G2;G3];
end
function xVec = state2FreeVariableVector(x, dx, ddx, h, parms)

xVec = zeros(parms.totalKnotNumber*parms.nVarSeg+2, 1);
for i = 1:parms.totalKnotNumber
    xSegment = [x(:, i); ...
        dx(:, i); ...
        ddx(:, i); ...
        ];
    xVec((i - 1)*parms.nVarSeg+(1:parms.nVarSeg), 1) = xSegment;
    
end

xVec(end-1) = h(1);
xVec(end) = h(2);
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

h(1) = aVec(end);
h(2) = aVec(end-1);
end

function xVec = initialGuess(parms)
% h
h = [3, 2];
x1 = [linspace(1, 1, parms.phase(1).knotNumber); ...
    linspace(parms.beta, parms.beta*2, parms.phase(1).knotNumber);];


% t2 = linspace(1,h(2),parms.phase(2).knotNumber);

x2 = [linspace(1, 1, parms.phase(2).knotNumber); ...
    linspace(parms.beta, parms.beta*2, parms.phase(2).knotNumber);];
x = [x1, x2];

% dx
dx = 0.1*ones(parms.ndof, parms.totalKnotNumber);
% ddx
ddx = 0.1*ones(parms.ndof, parms.totalKnotNumber);


xVec = state2FreeVariableVector(x, dx, ddx, h, parms);

end