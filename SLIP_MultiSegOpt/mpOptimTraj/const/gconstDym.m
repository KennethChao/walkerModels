function g = gconstDym(x,dx,ddx,parms)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


    iterRow = 1;
    iterCol = 1;
    oldInd = 0;
    
    for i=1:parms.phaseNum
        if i==1
            nRow = parms.ndof;
            nCol = parms.totalVarNumber;
            gI=nan( (nRow)*nCol,1);
            gJ=nan( (nRow)*nCol,1);
            gV=nan( (nRow)*nCol,1);    
        end
        
        if i ==1
            indexRange = 1: parms.phase(i).knotNumber;        
            xSeg = x(:,indexRange);
            dxSeg = dx(:,indexRange);      
        else        
            indexRange = parms.phase(i-1).knotNumber + (1: parms.phase(i).knotNumber);   
            xSeg = x(:,indexRange);
            dxSeg = dx(:,indexRange);                     
        end        
        
%         for j=1:2:(parms.phase(i).knotNumber)
        for j=1:(parms.phase(i).knotNumber)
                    
            gSegDym = parms.phase(i).jacobianDymFunc(xSeg(1,j),xSeg(2,j),dxSeg(1,j),dxSeg(2,j),parms.g, parms.k);

            sparseD = sparse(gSegDym);
            [SegI_D,SegJ_D,SegV_D] = find(sparseD);

            shiftInd = length(SegI_D);

            gI((1:shiftInd)+oldInd,1) = SegI_D+(nRow)*(iterRow-1);

%             gJ((1:shiftInd)+oldInd,1) = SegJ_D+(parms.nVarSeg*2)*(iterCol-1);
            gJ((1:shiftInd)+oldInd,1) = SegJ_D+(parms.nVarSeg)*(iterCol-1);

            gV((1:shiftInd)+oldInd,1) = SegV_D;

            oldInd = oldInd+shiftInd;
            iterRow = iterRow+1;
            iterCol = iterCol+1;
            if j == parms.phase(i).knotNumber
                iterCol = iterCol-0.5;
            end
        end
    end
    
%     g = sparse(gI(1:oldInd,1),gJ(1:oldInd,1),gV(1:oldInd,1),nRow*(parms.totalKnotNumber+2)/2,parms.totalVarNumber);
    g = sparse(gI(1:oldInd,1),gJ(1:oldInd,1),gV(1:oldInd,1),nRow*(parms.totalKnotNumber),parms.totalVarNumber);

end


