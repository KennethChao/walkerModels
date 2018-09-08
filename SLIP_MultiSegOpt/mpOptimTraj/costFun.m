function cost = costFun(x,dx,h,parms)
% cost = costFun(u)
%
cost = 0;
g = parms.g;
k = parms.k;
beta = parms.beta;

%% State in Cartesian space at the start of the stance phase
xStance0 = [x(:,1);dx(:,1)];
xStanceEnd = [x(:,parms.phase(1).knotNumber);dx(:,parms.phase(1).knotNumber)];

tFlight = getFlightTime(xStanceEnd(1),xStanceEnd(2),xStanceEnd(3),xStanceEnd(4),g,k,beta);

if isreal(tFlight) && tFlight>0

    cBoundary = boundaryCost(xStance0(1),xStance0(2),xStance0(3),xStance0(4),...
        xStanceEnd(1),xStanceEnd(2),xStanceEnd(3),xStanceEnd(4),....
        parms.g,parms.k,parms.beta );

else
    
    cBoundary = boundaryCostDummy(xStance0(1),xStance0(2),xStance0(3),xStance0(4),...
        xStanceEnd(1),xStanceEnd(2),xStanceEnd(3),xStanceEnd(4),....
        parms.g,parms.k,parms.beta );
    
end
%%
for i = 1:parms.phaseNum
    xSeg = x(:,(1:parms.phase(i).knotNumber)+parms.phase(i).x0knotNumber-1);
    dxSeg = dx(:,(1:parms.phase(i).knotNumber)+parms.phase(i).x0knotNumber-1);
    
%     if i == 1  
%         [~, xd, ~, zd] = polar2CartesianSLIP(xState( 1,:), dState( 1,:), xState( 2,:), dState( 2,:));
%         dCartesianState = [xd;zd];
%     else
%         dCartesianState = dState;
%     end
    cost = cost + sum(parms.phase(i).costFunc(xSeg(1,:),xSeg(2,:),dxSeg(1,:),dxSeg(2,:),h(i), g));
       
%     shiftIndex = shiftIndex + parms.phase(i).knotNumber;
end
cost = cost+parms.weightBoundary*sum(cBoundary);
end