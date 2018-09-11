function g = gconstPeriodicPattern(parms)
% cost = costFun(u)
%
% gcost = zeros(1,parms.totalVarNumber);
% g = parms.g;
shiftIndex = 0;

%         beta = parms.beta;
%         k = parms.k;

%%                   
%             parms.phase(i).jacobianBoundaryX0 = @jacobianBoundaryConstX0;
%             parms.phase(i).jacobianBoundaryXEnd = @jacobianBoundaryConstXEnd;
%     gcost(end-i) = parms.weightLagrangian*sum(parms.phase(i).jacobianCostH(xSeg(1,:),xSeg(2,:),dxSeg(1,:),dxSeg(2,:),h(i), g));

%     for j = (1:parms.phase(i).knotNumber )
%         gSeg = parms.weightLagrangian*parms.phase(i).jacobianCostX(xSeg(1,j),xSeg(2,j),dxSeg(1,j),dxSeg(2,j),h(i), g);
%         if j==1
%             gcost((1:parms.ndof*3) + shiftIndex) = gSeg + parms.weightBoundary*sum(gSegBoundary0);
%         elseif j == parms.phase(i).knotNumber
%             gcost((1:parms.ndof*3) + shiftIndex) = gSeg + parms.weightBoundary*sum(gSegBoundaryEnd);
%         else
%             gcost((1:parms.ndof*3) + shiftIndex) = gSeg;
%         end
%         shiftIndex = shiftIndex + parms.ndof*3;
%     end
    

%     shiftIndex = shiftIndex + parms.phase(i).knotNumber;
% end
%     iter = 1;
    oldInd = 0;
    
%     for i=1:parms.phaseNum
%         if i==1
            nRow = parms.nBoundaryConst;
            nCol = parms.totalVarNumber;
            gI=nan( (nRow)*nCol,1);
            gJ=nan( (nRow)*nCol,1);
            gV=nan( (nRow)*nCol,1);    
%         end
% 
%         x0 = x(:,1);
%         dx0 = dx(:,1);
% 
%         xEnd = x(:,end);
%         dxEnd = dx(:,end);  
%             
%   
% tFlight = getFlightTime(xEnd(1),xEnd(2),dxEnd(1),dxEnd(2),g,k,beta);


for i = 1:parms.phaseNum
            gSegBoundary0 = parms.phase(i).jacobianBoundaryCostX0Pattern();
            gSegBoundaryEnd = parms.phase(i).jacobianBoundaryCostXEndPattern();
    


%%                   
%             parms.phase(i).jacobianBoundaryX0 = @jacobianBoundaryConstX0;
%             parms.phase(i).jacobianBoundaryXEnd = @jacobianBoundaryConstXEnd;
%         x0 = xSeg(:,1);
%         dx0 = dxSeg(:,1);



            sparseB = sparse(gSegBoundary0);
            [SegI_B,SegJ_B,SegV_B] = find(sparseB);

            shiftInd = length(SegI_B);

            gI((1:shiftInd)+oldInd,1) = SegI_B;

            gJ((1:shiftInd)+oldInd,1) = SegJ_B+(parms.nVarSeg)*(parms.phase(i).x0knotNumber-1);                
            
            gV((1:shiftInd)+oldInd,1) = SegV_B;

            oldInd = oldInd+shiftInd;

 %%         

            sparseB = sparse(gSegBoundaryEnd);
            [SegI_B,SegJ_B,SegV_B] = find(sparseB);

            shiftInd = length(SegI_B);

            gI((1:shiftInd)+oldInd,1) = SegI_B;

            gJ((1:shiftInd)+oldInd,1) = SegJ_B+(parms.nVarSeg)*(parms.phase(i).x0knotNumber-1 + parms.phase(i).knotNumber-1);                
            
            gV((1:shiftInd)+oldInd,1) = SegV_B;
            
            oldInd = oldInd+shiftInd;      
            
            shiftInd = 2;
            
            gI((1:shiftInd)+oldInd,1) = [1;2];
            gJ((1:shiftInd)+oldInd,1) = [parms.totalVarNumber;parms.totalVarNumber];
            gV((1:shiftInd)+oldInd,1) = [1;1];
            
            oldInd = oldInd+shiftInd;   
    
    g = sparse(gI(1:oldInd,1),gJ(1:oldInd,1),gV(1:oldInd,1),parms.nBoundaryConst,parms.totalVarNumber);

% end

end