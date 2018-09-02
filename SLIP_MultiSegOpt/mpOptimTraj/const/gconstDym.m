function g = gconstDym(x,dx,ddx,parms)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


    iter = 1;
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
        
        for j=1:(parms.phase(i).knotNumber)
                    
            gSegDym = parms.phase(i).jacobianDymFunc(xSeg(1,j),xSeg(2,j),dxSeg(1,j),dxSeg(2,j),parms.g, parms.k);

            sparseD = sparse(gSegDym);
            [SegI_D,SegJ_D,SegV_D] = find(sparseD);

            shiftInd = length(SegI_D);

            gI((1:shiftInd)+oldInd,1) = SegI_D+(nRow)*(iter-1);

            gJ((1:shiftInd)+oldInd,1) = SegJ_D+(parms.nVarSeg)*(iter-1);

            gV((1:shiftInd)+oldInd,1) = SegV_D;

            oldInd = oldInd+shiftInd;
            iter = iter+1;
        
        end
    end
    
    g = sparse(gI(1:oldInd,1),gJ(1:oldInd,1),gV(1:oldInd,1),nRow*parms.totalKnotNumber,parms.totalVarNumber);

end


