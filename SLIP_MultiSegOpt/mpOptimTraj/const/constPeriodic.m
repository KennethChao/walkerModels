function c = constPeriodic(x, dx, sigma, parms)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
g = parms.g;
k = parms.k;
beta = parms.beta;

%% State in Cartesian space at the start of the stance phase
% xStance0 = [x(:,1);dx(:,1)];
xStance0 = [x(:,1);dx(:,1)];
xStanceEnd = [x(:,parms.phase(1).knotNumber);dx(:,parms.phase(1).knotNumber)];

tFlight = getFlightTime(xStanceEnd(1),xStanceEnd(2),xStanceEnd(3),xStanceEnd(4), xStance0(2),g);

if isreal(tFlight) && tFlight>0
%     tFlight
    cBoundary = boundaryCost(xStance0(1),xStance0(2),xStance0(3),xStance0(4),...
        xStanceEnd(1),xStanceEnd(2),xStanceEnd(3),xStanceEnd(4),....
        parms.g,parms.k,parms.beta );

else
    cBoundary = boundaryCostDummy(xStance0(1),xStance0(2),xStance0(3),xStance0(4),...
        xStanceEnd(1),xStanceEnd(2),xStanceEnd(3),xStanceEnd(4),....
        parms.g,parms.k,parms.beta );
    
end

c = cBoundary-sigma;
end
