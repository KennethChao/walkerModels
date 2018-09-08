function g = gconstBoundary(x,dx,ddx,parms)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


    iter = 1;
    oldInd = 0;
    
    for i=1:parms.phaseNum
        if i==1
            nRow = parms.nBoundaryConst;
            nCol = parms.totalVarNumber;
            gI=nan( (nRow)*nCol,1);
            gJ=nan( (nRow)*nCol,1);
            gV=nan( (nRow)*nCol,1);    
        end
        
      
        indexRange = parms.phase(i).x0knotNumber-1 + (1: parms.phase(i).knotNumber);   
        xSeg = x(:,indexRange);
        dxSeg = dx(:,indexRange);                     
                
        beta = parms.beta;
        g = parms.g;
        k = parms.k;

%%                   
%             parms.phase(i).jacobianBoundaryX0 = @jacobianBoundaryConstX0;
%             parms.phase(i).jacobianBoundaryXEnd = @jacobianBoundaryConstXEnd;
        x0 = xSeg(:,1);
        dx0 = dxSeg(:,1);



            gSegBoundary0 = parms.phase(i).jacobianBoundaryConstX0(...
                x0(1),x0(2),dx0(1),dx0(2), g, k, beta);


            sparseB = sparse(gSegBoundary0);
            [SegI_B,SegJ_B,SegV_B] = find(sparseB);

            shiftInd = length(SegI_B);

            gI((1:shiftInd)+oldInd,1) = SegI_B;

            gJ((1:shiftInd)+oldInd,1) = SegJ_B+(parms.nVarSeg)*(parms.phase(i).x0knotNumber-1);                
            
            gV((1:shiftInd)+oldInd,1) = SegV_B;

            oldInd = oldInd+shiftInd;

 %%         
%  if isreal(tFlight) && tFlight>0             
%             gSegBoundaryEnd = parms.phase(i).jacobianBoundaryXEnd(...
%                 xEnd(1),xEnd(2),dxEnd(1),dxEnd(2), g,k, beta);
%  else
%             gSegBoundaryEnd = parms.phase(i).jacobianBoundaryXEndFake(...
%                 xEnd(1),xEnd(2),dxEnd(1),dxEnd(2), g,k, beta);
%  end
%             sparseB = sparse(gSegBoundaryEnd);
%             [SegI_B,SegJ_B,SegV_B] = find(sparseB);
% 
%             shiftInd = length(SegI_B);
% 
%             gI((1:shiftInd)+oldInd,1) = SegI_B;
% 
%             gJ((1:shiftInd)+oldInd,1) = SegJ_B+(parms.nVarSeg)*(parms.phase(i).x0knotNumber-1 + parms.phase(i).knotNumber-1);                
%             
%             gV((1:shiftInd)+oldInd,1) = SegV_B;
%             
%             oldInd = oldInd+shiftInd;                
    

  
            
            
 
    end
    
    g = sparse(gI(1:oldInd,1),gJ(1:oldInd,1),gV(1:oldInd,1),parms.nBoundaryConst,parms.totalVarNumber);

end


