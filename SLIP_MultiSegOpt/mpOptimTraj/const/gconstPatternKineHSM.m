function g = gconstPatternKineHSM(parms) %x, dx, ddx, h,
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    ndof = parms.ndof;
    
        gSeg1 = [  eye(ndof), eye(ndof), zeros(ndof),...
                 zeros(ndof), eye(ndof), zeros(ndof),...
                   eye(ndof), eye(ndof), zeros(ndof)];

        gSeg2 = [zeros(ndof),    eye(ndof),  eye(ndof),... 
                 zeros(ndof),  zeros(ndof),  eye(ndof),... 
                 zeros(ndof),    eye(ndof),  eye(ndof)];               

        gSeg3 = [ eye(ndof),   eye(ndof), zeros(ndof),...
                  eye(ndof), zeros(ndof), zeros(ndof),...
                  eye(ndof),   eye(ndof), zeros(ndof)];

        gSeg4 = [zeros(ndof),  eye(ndof), eye(ndof),... 
                 zeros(ndof),  eye(ndof), zeros(ndof),... 
                 zeros(ndof),  eye(ndof), eye(ndof),];   

        gSeg1h = ones(ndof,1);

        gSeg2h = ones(ndof,1);

        gSeg3h = ones(ndof,1);             

        gSeg4h = ones(ndof,1);                          

        gSegKine = [gSeg1;
                    gSeg2;
                    gSeg3;
                    gSeg4];

        gKineh =  [gSeg1h;
                   gSeg2h;
                   gSeg3h;
                   gSeg4h];    
    iter = 1;
    oldInd = 0;
    shiftIndex = 0;
    for i=1:parms.phaseNum
        if i==1
            nRow = parms.ndof*4;
            gI=nan( (nRow)*(parms.totaHSMCnstNumber),1);
            gJ=nan( (nRow)*(parms.totaHSMCnstNumber),1);
            gV=nan( (nRow)*(parms.totaHSMCnstNumber),1);    
        end
        
        for j=1:2:(parms.phase(i).knotNumber-2)
            
        sparseK = sparse(gSegKine);
        [SegI_K,SegJ_K,SegV_K] = find(sparseK);
              
        sparseKh = sparse(gKineh);
        [SegI_Kh,~,SegV_Kh] = find(sparseKh);            

        shiftInd = size(SegI_K,1)+size(SegI_Kh,1);

        gI((1:shiftInd)+oldInd,1) = [SegI_K+(nRow)*(iter-1);...
                                   SegI_Kh+(nRow)*(iter-1)];

        gJ((1:shiftInd)+oldInd,1) = [SegJ_K + (parms.nVarSeg)*((iter-1 )*2 + (i-1));...
                                       (parms.totalVarNumber-i+1)*ones((size(SegI_Kh,1)),1)];

        gV((1:shiftInd)+oldInd,1) = [SegV_K;...
                                   SegV_Kh];
        oldInd = oldInd+shiftInd;
        iter = iter+1;
        
        end
        
        shiftIndex = shiftIndex + parms.phase(i).knotNumber;
    end
    
    g = sparse(gI(1:oldInd,1),gJ(1:oldInd,1),gV(1:oldInd,1),nRow*parms.totaHSMCnstNumber,parms.totalVarNumber);

end

