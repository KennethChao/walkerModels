function g = gconstKineHSM( parms) %x, dx, ddx, h,
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    x1 = [linspace(1, 1, parms.phase(1).knotNumber); ...
    linspace(parms.beta, parms.beta*2, parms.phase(1).knotNumber);];


    % t2 = linspace(1,h(2),parms.phase(2).knotNumber);

    x2 = [linspace(1, 1, parms.phase(2).knotNumber); ...
        linspace(parms.beta, parms.beta*2, parms.phase(2).knotNumber);];

    x = [x1, x2];

    % dx
    dx = ones(parms.ndof, parms.totalKnotNumber);
    % ddx
    ddx = zeros(parms.ndof, parms.totalKnotNumber);

    h  = [1;2];

    iter = 1;

    oldInd = 0;

    for i=1:length(h)
        if i==1
              nRow = parms.ndof*4;
              nCol = parms.nVarSeg*3;
              gI=nan( (nRow)*(nCol+1)*(parms.totalKnotNumber-4),1);
              gJ=nan( (nRow)*(nCol+1)*(parms.totalKnotNumber-4),1);
              gV=nan( (nRow)*(nCol+1)*(parms.totalKnotNumber-4),1);    
        end        
        
        for j=1:(parms.phase(i).knotNumber-2)
                                                 
        gSeg1 = [-1*eye(2), -1/6*h(i)*eye(2), zeros(2),...
                  zeros(2), -4/6*h(i)*eye(2), zeros(2),...
                    eye(2),  1/6*h(i)*eye(2), zeros(2)];
               
        gSeg2 = [zeros(2),  -1*eye(2), -1/6*h(i)*eye(2),... 
                 zeros(2),   zeros(2), -4/6*h(i)*eye(2),... 
                 zeros(2),     eye(2),  1/6*h(i)*eye(2)];               

        gSeg3 = [-1/2*eye(2), -1/8*h(i)*eye(2), zeros(2),...
                      eye(2),         zeros(2), zeros(2),...
                 -1/2*eye(2),  1/8*h(i)*eye(2), zeros(2)];
               
        gSeg4 = [zeros(2),  -1/2*eye(2), -1/8*h(i)*eye(2),... 
                 zeros(2),       eye(2),         zeros(2),... 
                 zeros(2),  -1/2*eye(2),  1/8*h(i)*eye(2),];   
             
        gSeg1h = -1/6*dx(:,j)-4/6*dx(:,j+1)-1/6*dx(:,j+2);
             
        gSeg2h = -1/6*ddx(:,j)-4/6*ddx(:,j+1)-1/6*ddx(:,j+2);
             
        gSeg3h = -1/8*dx(:,j)+ zeros(2,1)+1/8*dx(:,j+2);             
             

        gSeg4h = -1/8*ddx(:,j)+zeros(2,1)+1/8*ddx(:,j+2);                          

     

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

        gJ((1:shiftInd)+oldInd,1) = [SegJ_K+(parms.nVarSeg)*(iter-1+(i-1)*2);...
                                       (parms.totalVarNumber-i+1)*ones((size(SegI_Kh,1)),1)];

        gV((1:shiftInd)+oldInd,1) = [SegV_K;...
                                   SegV_Kh];
        oldInd = oldInd+shiftInd;
        iter = iter+1;
%             if j==(parms.phase(i).knotNumber-2) && i == length(h)-1
%                 iter = iter+1;
%             end
        end
    end
    
    g = sparse(gI(1:oldInd,1),gJ(1:oldInd,1),gV(1:oldInd,1),nRow*(parms.totalKnotNumber-4),parms.totalVarNumber);

end
