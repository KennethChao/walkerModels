function gcost = gcostFun(x,dx,h,parms)
% cost = costFun(u)
%
gcost = zeros(1,parms.totalVarNumber);
shiftIndex = 0;
for i = 1:parms.phaseNum
    xSeg = x(:,(1:parms.phase(i).knotNumber)+parms.phase(i).x0knotNumber-1);
    dxSeg = dx(:,(1:parms.phase(i).knotNumber)+parms.phase(i).x0knotNumber-1);
    
%     if i == 1  
%         [~, xd, ~, zd] = polar2CartesianSLIP(xState( 1,:), dState( 1,:), xState( 2,:), dState( 2,:));
%         dCartesianState = [xd;zd];
%     else
%         dCartesianState = dState;
%     end
    gcost(end-i+1) = sum(parms.phase(i).jacobianCostH(xSeg(1,:),xSeg(2,:),dxSeg(1,:),dxSeg(2,:),h(i)));

    for j = (1:parms.phase(i).knotNumber )
        gcost((1:parms.ndof*3) + shiftIndex) = parms.phase(i).jacobianCostX(xSeg(1,j),xSeg(2,j),dxSeg(1,j),dxSeg(2,j),h(i));
        shiftIndex = shiftIndex + parms.ndof*3;
    end
    
    
%     shiftIndex = shiftIndex + parms.phase(i).knotNumber;
end

end