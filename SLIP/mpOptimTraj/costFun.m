function cost = costFun(x,dx,h,parms)
% cost = costFun(u)
%
cost = 0;
shiftIndex = 0;

for i = 1:phaseNum
    xState = x(:,(1:parms.phase(i).knotNumber)+shiftIndex);
    dState = dx(:,(1:parms.phase(i).knotNumber)+shiftIndex);
    
    if i == 1  
        [~, xd, ~, zd] = polar2CartesianSLIP(xState( 1,:), dState( 1,:), xState( 2,:), dState( 2,:));
        dCartesianState = [xd;zd];
    else
        dCartesianState = dState;
    end
    cost = cost + dCartesianState.^2*h(i)*parms.phase(i).knotNumber;
    
    shiftIndex = shiftIndex + parms.phase(i).knotNumber;
end

end