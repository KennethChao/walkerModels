function gcost = gcostFun(x,dx,h,parms)
% cost = costFun(u)
%
gcost = zeros(1,parms.totalVarNumber);
g = parms.g;
shiftIndex = 0;

        beta = parms.beta;
        k = parms.k;

%%                   
%             parms.phase(i).jacobianBoundaryX0 = @jacobianBoundaryConstX0;
%             parms.phase(i).jacobianBoundaryXEnd = @jacobianBoundaryConstXEnd;
        x0 = x(:,1);
        dx0 = dx(:,1);

        xEnd = x(:,end);
        dxEnd = dx(:,end);  
            
  
tFlight = getFlightTime(xEnd(1),xEnd(2),dxEnd(1),dxEnd(2),g,k,beta);


for i = 1:parms.phaseNum
    xSeg = x(:,(1:parms.phase(i).knotNumber)+parms.phase(i).x0knotNumber-1);
    dxSeg = dx(:,(1:parms.phase(i).knotNumber)+parms.phase(i).x0knotNumber-1);
    
    if isreal(tFlight) && tFlight>0  
            gSegBoundary0 = parms.phase(i).jacobianBoundaryCostX0(...
                x0(1),x0(2),dx0(1),dx0(2), xEnd(1),xEnd(2),dxEnd(1),dxEnd(2), g, k, beta);
            
            gSegBoundaryEnd = parms.phase(i).jacobianBoundaryCostXEnd(...
                x0(1),x0(2),dx0(1),dx0(2), xEnd(1),xEnd(2),dxEnd(1),dxEnd(2), g,k, beta);
    else
            gSegBoundary0 = parms.phase(i).jacobianBoundaryCostX0Dummy(...
                x0(1),x0(2),dx0(1),dx0(2), xEnd(1),xEnd(2),dxEnd(1),dxEnd(2), g, k, beta);    

            gSegBoundaryEnd = parms.phase(i).jacobianBoundaryCostXEndDummy(...
                x0(1),x0(2),dx0(1),dx0(2), xEnd(1),xEnd(2),dxEnd(1),dxEnd(2), g,k, beta);            
            
    end
    
    gcost(end-i) = parms.weightLagrangian*sum(parms.phase(i).jacobianCostH(xSeg(1,:),xSeg(2,:),dxSeg(1,:),dxSeg(2,:),h(i), g));

    for j = (1:parms.phase(i).knotNumber )
        gSeg = parms.weightLagrangian*parms.phase(i).jacobianCostX(xSeg(1,j),xSeg(2,j),dxSeg(1,j),dxSeg(2,j),h(i), g);
        if j==1
            gcost((1:parms.ndof*3) + shiftIndex) = gSeg + parms.weightBoundary*sum(gSegBoundary0);
        elseif j == parms.phase(i).knotNumber
            gcost((1:parms.ndof*3) + shiftIndex) = gSeg + parms.weightBoundary*sum(gSegBoundaryEnd);
        else
            gcost((1:parms.ndof*3) + shiftIndex) = gSeg;
        end
        shiftIndex = shiftIndex + parms.ndof*3;
    end
    

%     shiftIndex = shiftIndex + parms.phase(i).knotNumber;
end

end