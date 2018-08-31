function g = gconstDymPattern(parms)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


    iter = 1;
    oldInd = 0;
    
    for i=1:parms.phaseNum
        if i==1
            nRow = parms.ndof;
            gI=nan( (nRow)*(parms.totalKnotNumber),1);
            gJ=nan( (nRow)*(parms.totalKnotNumber),1);
            gV=nan( (nRow)*(parms.totalKnotNumber),1);    
        end
         
        
        for j=1:(parms.phase(i).knotNumber)
                    
%             gSegDym = parms.phase(i).dymFunc(xSeg(1,i),xSeg(2,i),dxSeg(1,i),dxSeg(2,i),parms.g, parms.k);
            gSegDym = [ones(2,4),eye(2)];
            gSegDym = [ones(2,6)];
            sparseD = sparse(gSegDym);
            [SegI_D,SegJ_D,SegV_D] = find(sparseD);

            shiftInd = size(SegI_D,1);

            gI((1:shiftInd)+oldInd,1) = SegI_D+(nRow)*(iter-1);

            gJ((1:shiftInd)+oldInd,1) = SegJ_D+(parms.nVarSeg)*(iter-1 );

            gV((1:shiftInd)+oldInd,1) = SegV_D;

            oldInd = oldInd+shiftInd;
            iter = iter+1;
        
        end
    end
    
    g = sparse(gI(1:oldInd,1),gJ(1:oldInd,1),gV(1:oldInd,1),nRow*parms.totalKnotNumber,parms.totalVarNumber);

end


