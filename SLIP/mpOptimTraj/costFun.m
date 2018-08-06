function cost = costFun(x,dx,h,parms)
% cost = costFun(u)
%

cost = 0;
shiftIndex = 0;

for i = 1:phaseNum
    xState = x(:,(1:parms.phase(i).knotNumber)+shiftIndex);
    dState = dx(:,(1:parms.phase(i).knotNumber)+shiftIndex);

    
    if i == 1  
        dCartesianState = ...
        [-dState(1,: ) .* cos(xState(2,: )) + xState( 1,:).* dState(2,:) .* sin(xState(2,:)); ...
         dState(1,: ) .* sin(xState(2,: )) + xState( 1,:).* dState(2,:) .* cos(xState(2,: ))];
        
    else
       dCartesianState = dState;
    end
    cost = cost + dCartesianState.^2*h(i)*parms.phase(i).knotNumber;
    
    shiftIndex = shiftIndex + parms.phase(i).knotNumber;
end



end