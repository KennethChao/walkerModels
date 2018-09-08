function g = gconstBoundaryPattern(parms)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


%     iter = 1;
    oldInd = 0;
    
    for i=1:parms.phaseNum
        if i==1
            nRow = parms.nBoundaryConst;
            nCol = parms.totalVarNumber;
            gI=nan( (nRow)*nCol,1);
            gJ=nan( (nRow)*nCol,1);
            gV=nan( (nRow)*nCol,1);    
        end
        
%%                   
%             parms.phase(i).jacobianBoundaryX0 = @jacobianBoundaryConstX0;
%             parms.phase(i).jacobianBoundaryXEnd = @jacobianBoundaryConstXEnd;

            gSegBoundary0 = parms.phase(i).jacobianBoundaryConstX0Pattern();

            sparseB = sparse(gSegBoundary0);
            [SegI_B,SegJ_B,SegV_B] = find(sparseB);

            shiftInd = length(SegI_B);

            gI((1:shiftInd)+oldInd,1) = SegI_B;

            gJ((1:shiftInd)+oldInd,1) = SegJ_B+(parms.nVarSeg)*(parms.phase(i).x0knotNumber-1);                
            
            gV((1:shiftInd)+oldInd,1) = SegV_B;

            oldInd = oldInd+shiftInd;

%%          
%             gSegBoundaryEnd = parms.phase(i).jacobianBoundaryXEndPattern();
% 
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


