function g = gconstKineHSM(x, dx, ddx, h, parms) %x, dx, ddx, h,
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
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
        
        if i ==1
            indexRange = 1: parms.phase(i).knotNumber;        
            dxSeg = dx(:,indexRange);
            ddxSeg = ddx(:,indexRange);           
        else        
            indexRange = parms.phase(i-1).knotNumber + (1: parms.phase(i).knotNumber);   
            dxSeg = dx(:,indexRange);
            ddxSeg = ddx(:,indexRange);                      
        end        
        
        for j=1:2:(parms.phase(i).knotNumber-2)
            
          
                                                 
        gSeg1 = [-1*eye(2), -1/6*h(i)*eye(2), zeros(2),...
                  zeros(2), -4/6*h(i)*eye(2), zeros(2),...
                    eye(2), -1/6*h(i)*eye(2), zeros(2)];
               
        gSeg2 = [zeros(2),  -1*eye(2), -1/6*h(i)*eye(2),... 
                 zeros(2),   zeros(2), -4/6*h(i)*eye(2),... 
                 zeros(2),     eye(2), -1/6*h(i)*eye(2)];               

        gSeg3 = [-1/2*eye(2), -1/8*h(i)*eye(2), zeros(2),...
                      eye(2),         zeros(2), zeros(2),...
                 -1/2*eye(2),  1/8*h(i)*eye(2), zeros(2)];
               
        gSeg4 = [zeros(2),  -1/2*eye(2), -1/8*h(i)*eye(2),... 
                 zeros(2),       eye(2),         zeros(2),... 
                 zeros(2),  -1/2*eye(2),  1/8*h(i)*eye(2),];   
             
        gSeg1h = -1/6*dxSeg(:,j)-4/6*dxSeg(:,j+1)-1/6*dxSeg(:,j+2);
             
        gSeg2h = -1/6*ddxSeg(:,j)-4/6*ddxSeg(:,j+1)-1/6*ddxSeg(:,j+2);
             
        gSeg3h = -1/8*dxSeg(:,j)+1/8*dxSeg(:,j+2);             
             

        gSeg4h = -1/8*ddxSeg(:,j)+1/8*ddxSeg(:,j+2);                          

        
        gSegKine = [gSeg1;
                    gSeg2;
                    gSeg3;
                    gSeg4];

        gKineh =  [gSeg1h;
                   gSeg2h;
                   gSeg3h;
                   gSeg4h];

        sparseK = sparse(gSegKine);
        [SegI_K,SegJ_K,SegV_K] = find(sparseK);
              
        sparseKh = sparse(gKineh);
        [SegI_Kh,~,SegV_Kh] = find(sparseKh);        
        
        shiftInd = size(SegI_K,1)+size(SegI_Kh,1);

        gI((1:shiftInd)+oldInd,1) = [SegI_K+(nRow)*(iter-1);...
                                   SegI_Kh+(nRow)*(iter-1)];

        gJ((1:shiftInd)+oldInd,1) = [SegJ_K+(parms.nVarSeg)*((iter-1 )*2 + (i-1));...
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

