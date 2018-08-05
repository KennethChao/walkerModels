%MAIN.m  --  simple walker trajectory optimization
%
% This script sets up a trajectory optimization problem for a simple model
% of walking, and solves it using OptimTraj. The walking model is a double
% pendulum, with point feet, no ankle torques, impulsive heel-strike (but
% not push-off), and continuous hip torque. Both legs have inertia. Cost
% function is minimize integral of torque-squared.
%
%

clc; clear;
% addpath ../../

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                  Parameters for the dynamics function                   %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%% Parameter Set

parms.g = 0.05;
parms.beta = 72/180*pi;
parms.k = 12;

parms.delta0 = 0.1;

%
parms.phase(1).knotNumber = 20;
parms.phase(1).dymFunc = @dymModelStanceDimensionless;

parms.phase(2).knotNumber = 10;

totalKnotNumber = 0;
for i = length(parms.phase)
   totalKnotNumber = totalKnotNumber + parms.phase(i).knotNumber; 
end
parms.totalKnotNumber = totalKnotNumber;
qq=0;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                       Set up function handles                           %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.func.dynamics = @(t,x,u)( dymModelStanceDimensionless(t, x,u,parms) );

problem.func.pathObj = @(t,x,u)( cost(t,x) );

problem.func.bndCst = @(t0,x0,tF,xF)( periodicGait(xF,x0,parms) );


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%               Set up bounds on time, state, and control                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
t0 = 0;  tF = 0.5;
problem.bounds.initialTime.low = t0;
problem.bounds.initialTime.upp = t0;
problem.bounds.finalTime.low = tF;
problem.bounds.finalTime.upp = tF+1;

% State: [q1;q2;dq1;dq2];

problem.bounds.state.low = [0.01 ;-inf;parms.beta-0.1; -inf];
problem.bounds.state.upp = [1.01; inf;  pi;  inf];

% stepAngle = 0.2;
problem.bounds.initialState.low = [1;-inf; parms.beta; -inf];
problem.bounds.initialState.upp = [1;inf; parms.beta;  inf];

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%              Create an initial guess for the trajectory                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% For now, just assume a linear trajectory between boundary values

problem.guess.time = [t0, tF];

% stepRate = (2*stepAngle)/(tF-t0);
x0 = [1; -1;  parms.beta; 0.1];
xF = [1; 1;  parms.beta*2; 0.1];
problem.guess.state = [x0, xF];

problem.guess.control = [0, 0];


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           Options:                                      %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%


%NOTE:  Here I choose to run the optimization twice, mostly to demonstrate
%   functionality, although this can be important on harder problems. I've
%   explicitly written out many options below, but the solver will fill in
%   almost all defaults for you if they are ommitted.

% method = 'trapezoid';
method = 'hermiteSimpson';
% method = 'chebyshev';
% method = 'rungeKutta';
% method = 'gpops';

switch method
    case 'trapezoid'
        
        % First iteration: get a more reasonable guess
        problem.options(1).nlpOpt = optimset(...
            'Display','iter',...   % {'iter','final','off'}
            'TolFun',1e-3,...
            'MaxFunEvals',1e4);   %options for fmincon
        problem.options(1).verbose = 3; % How much to print out?
        problem.options(1).method = 'trapezoid'; % Select the transcription method
        problem.options(1).trapezoid.nGrid = 10;  %method-specific options
        
        
        % Second iteration: refine guess to get precise soln
        problem.options(2).nlpOpt = optimset(...
            'Display','iter',...   % {'iter','final','off'}
            'TolFun',1e-6,...
            'MaxFunEvals',5e4);   %options for fmincon
        problem.options(2).verbose = 3; % How much to print out?
        problem.options(2).method = 'trapezoid'; % Select the transcription method
        problem.options(2).trapezoid.nGrid = 25;  %method-specific options
        
    case 'hermiteSimpson'
        
        % First iteration: get a more reasonable guess
        problem.options(1).nlpOpt = optimset(...
            'Display','iter',...   % {'iter','final','off'}
            'TolFun',1e-3,...
            'MaxFunEvals',1e4);   %options for fmincon
        problem.options(1).verbose = 3; % How much to print out?
        problem.options(1).method = 'hermiteSimpson'; % Select the transcription method
        problem.options(1).hermiteSimpson.nSegment = 6;  %method-specific options
        
        
        % Second iteration: refine guess to get precise soln
        problem.options(2).nlpOpt = optimset(...
            'Display','iter',...   % {'iter','final','off'}
            'TolFun',1e-6,...
            'MaxFunEvals',5e4);   %options for fmincon
        problem.options(2).verbose = 3; % How much to print out?
        problem.options(2).method = 'hermiteSimpson'; % Select the transcription method
        problem.options(2).hermiteSimpson.nSegment = 15;  %method-specific options
        

end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           Solve!                                        %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%%%%% THE KEY LINE:
soln = optimTraj(problem);
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
function c = cost(xVec, parms)

%%
% t = 
% x = 
% dx = 
% ddx =
%%
c = costFun(t,dx);
end %function end

function c = constAll(xVec, parms)

%%
% x = 
% dx = 
% ddx =
%%
c0=constKineHSM(x,dx,ddx, parms);
c1=constDymHSM(aVec, aInd, parms);
c2=constPeriodic(aVec, aInd, parms);
c = [c0;c1;c2];
end %function end

% function g = gconstAll(aVec, aInd,cfg)
% g0=gKineHSM(aVec, aInd, cfg);
% g1=gDymHSM(aVec, aInd, cfg);
% g2=gContact(aVec, aInd, cfg);
% g3=gPeriodic(aVec, aInd, cfg);
% g = [g0;g1;g2;g3];%g1;g2;g3
% % g = [g0];%g1;g2;g3
% end %function end
% 
% function Pattern = GP(cfg)  
% G0 = gKineHSMPat(cfg);
% G1 = gDymHSMPat(cfg);
% G2 = gContactPat(cfg);
% G3 = gPeriodicPat(cfg);
% Pattern = [G0;G1;G2;G3];
% end
