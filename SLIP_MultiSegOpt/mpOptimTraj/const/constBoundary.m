function c = constBoundary(x, dx, parms)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
g = parms.g;
k = parms.k;
beta = parms.beta;

%% State in Cartesian space at the start of the stance phase
xStance0 = [x(:,1);dx(:,1)];


    c = boundaryConst(xStance0(1),xStance0(2),xStance0(3),xStance0(4),...
        parms.g,parms.k,parms.beta );



end
