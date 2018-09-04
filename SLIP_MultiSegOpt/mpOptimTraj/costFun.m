function cost = costFun(x,dx,h,parms)
% cost = costFun(u)
%
cost = 0;
g = parms.g;

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

end