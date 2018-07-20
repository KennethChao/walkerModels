%%
close all;


%% Initial Condition (Guessed)
k = 15;
delta = 0.1392;

betaVec = 66:2:74
for k = 1:length(betaVec)

beta = betaVec(k)/180*pi;
% beta = 66/180*pi;
g = 0.46;

%% Parameter Set

parms = {};
parms.g = g;
parms.beta = beta;
parms.k = k;
parms.delta = delta;

options = optimset('Display','off');
% options = optimset('Display','iter');

%% Sampling Number and Result Buffer
samplingNumbK = 25;
samplingNumbDelta = 10;

stableSolution = nan*zeros(samplingNumbDelta,samplingNumbK);
unstableSolution = nan*zeros(samplingNumbDelta,samplingNumbK);

kVec = linspace(0.1,25,samplingNumbK);
deltaVec = linspace(1e-5,1.5,samplingNumbDelta);

%%

for i = 1:samplingNumbK
    parms.k = kVec(i);
    for j = 1:samplingNumbDelta
    delta0 = deltaVec(j);
           
    [x,fval,exitflag,output]  =fminunc(@(x)oneStepSimulation(x,parms),delta0,options);
        if  fval<1e-2 && x>1e-5
            
        eigenValue =  perturbedEigenValue(x,parms);
        
            if abs(eigenValue)<1
                stableSolution(j,i) = x;
            else
                unstableSolution(j,i) = x;
            end
        end
    end
        msg=sprintf('%d th iteration for k',i);
        disp(msg)
end


% plot(solution','bo')
for i = 1:3
    plot(kVec,stableSolution(i,:),'bo')
    plot(kVec,unstableSolution(i,:),'ro')
    hold on
end

axis([0 25, 0 1])
end
%%
function result =  oneStepSimulation(delta0,parms)
beta = parms.beta;
% delta = parms.delta;
% perturbedDelta = delta+perturbationl;
deltaVec = delta0;
delta =delta0;
normVel = 1;
for i=1:1
dymStance = @(t,x) dymModelStance(t,x, parms);

tspan = 0:0.01:100;
options = odeset('Event',@slipEventFcn);
x0 = [1,-normVel*cos(beta-delta),beta,normVel*sin(beta-delta)];
[t,x,te,xe,ie] = ode45(dymStance,tspan,x0,options);     % Runge-Kutta 4th/5th order ODE solver
% te
% figure()
% plot(x(:,[1 3]))
xVec = x(:,1).*cos(x(:,3))*-1;
yVec = x(:,1).*sin(x(:,3));

% figure()
% plot(xVec,yVec);

y0 = x(end,1)*sin(x(end,3));
yd0 = x(end,2)*cos(x(end,3)-pi/2) - x(end,1)*x(end,4)*sin(x(end,3)-pi/2);
x0 = [y0,yd0];
options = odeset('Event',@(t,x)tdEventFcn(t,x,beta));
dymFlight = @(t,a) dymModelFlight(t,a, parms);
[t,x2,te2,ae2,ie2] = ode45(dymFlight,tspan,x0,options);     % Runge-Kutta 4th/5th order ODE solver
% te2
% figure()
% plot(x2(:,[1]))
velVec = [-x2(end,2),x(end,2)*sin(x(end,3)-pi/2) + x(end,1)*x(end,4)*cos(x(end,3)-pi/2)];
delta = atan2(velVec(1),velVec(2));
deltaVec = [deltaVec,delta];

normVel = norm(velVec);
end

diffDelta = deltaVec(end) -deltaVec(end-1) ;
diffBeta = 1-  norm(velVec);



result = abs(diffDelta) + abs(diffBeta);

% deltaVecShifted = deltaVec-deltaVec(end);
% eVar = deltaVecShifted(2:end)*pinv(deltaVecShifted(1:end-1));
% if eVar<1    
%     eVarBuffer = [eVar eVar];
% end
%     j
% plot(deltaVec)
end

function result =  perturbedEigenValue(delta0,parms)
beta = parms.beta;
% delta = parms.delta;
perturbation = 0.1;
perturbedDelta = delta0+perturbation;
delta =perturbedDelta;
normVel = 1;
for i=1:1
dymStance = @(t,x) dymModelStance(t,x, parms);

tspan = 0:0.01:100;
options = odeset('Event',@slipEventFcn);
x0 = [1,-normVel*cos(beta-delta),beta,normVel*sin(beta-delta)];
[t,x,te,xe,ie] = ode45(dymStance,tspan,x0,options);     % Runge-Kutta 4th/5th order ODE solver
% te
% figure()
% plot(x(:,[1 3]))
xVec = x(:,1).*cos(x(:,3))*-1;
yVec = x(:,1).*sin(x(:,3));

% figure()
% plot(xVec,yVec);

y0 = x(end,1)*sin(x(end,3));
yd0 = x(end,2)*cos(x(end,3)-pi/2) - x(end,1)*x(end,4)*sin(x(end,3)-pi/2);
x0 = [y0,yd0];
options = odeset('Event',@(t,x)tdEventFcn(t,x,beta));
dymFlight = @(t,a) dymModelFlight(t,a, parms);
[t,x2,te2,ae2,ie2] = ode45(dymFlight,tspan,x0,options);     % Runge-Kutta 4th/5th order ODE solver
% te2
% figure()
% plot(x2(:,[1]))
velVec = [-x2(end,2),x(end,2)*sin(x(end,3)-pi/2) + x(end,1)*x(end,4)*cos(x(end,3)-pi/2)];
delta = atan2(velVec(1),velVec(2));

normVel = norm(velVec);
end

result = (delta - perturbedDelta)/perturbation;

end


function dx = dymModelStance(t,x,parms)
%
g = parms.g;
k = parms.k;
%
l = x(1);
ld = x(2);
theta = x(3);
thetad = x(4);

ldd = l*thetad^2-k*(l-1)-g*sin(theta);
thetadd = (-g*l*cos(theta)-2*l*ld*thetad)/l^2;

dx = [ld;ldd;thetad;thetadd];

end

function dx = dymModelFlight(t,x,parms)
g = parms.g;
yd = x(2);
ydd = -g;

dx = [yd;ydd];

end


function [position,isterminal,direction] = slipEventFcn(t,x)
position = x(1)-1; % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 1;   % The zero can be approached from either direction
end

function [position,isterminal,direction] = tdEventFcn(t,x,beta)
position = x(1)-sin(beta); % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = -1;   % The zero can be approached from either direction
end