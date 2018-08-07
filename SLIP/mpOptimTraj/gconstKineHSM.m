function [outputArg1,outputArg2] = gconstKineHSM(inputArg1,parms)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    [n,m]=size(aVec);
    if n<m
        aVec=aVec'; 
    end
    h = aVec(end);
    
    iter = 1;

    oldInd = 0;

    for i=1:2:(cfg.nSeg-2)
        
        x1 = aVec(a.Xind(i,1):a.Xind(i,2));
        x2 = aVec(a.Xind(i+1,1):a.Xind(i+1,2));   
        x3 = aVec(a.Xind(i+2,1):a.Xind(i+2,2)); 
        xKine = [x1;x2;x3;h];        

        gSegKine = [PKine(xKine);...;
                    PKine2(xKine);...;
                    PKine3(xKine);...;
                    PKine4(xKine)];

        gKineh =  [PKineh(xKine);...
                   PKine2h(xKine);...
                   PKine3h(xKine);...
                   PKine4h(xKine)];
        if i==1
              nRow = cfg.ndof*4;
              nCol = cfg.nVarSeg*3;
              gI=nan( (nRow)*(nCol)*(cfg.nSeg-1)/2,1);
              gJ=nan( (nRow)*(nCol)*(cfg.nSeg-1)/2,1);
              gV=nan( (nRow)*(nCol)*(cfg.nSeg-1)/2,1);    
        end
              sparseK = sparse(gSegKine);
              [SegI_K,SegJ_K,SegV_K] = find(sparseK);
              
              sparseKh = sparse(gKineh);
              [SegI_Kh,~,SegV_Kh] = find(sparseKh);


              shiftInd = size(SegI_K,1)+size(SegI_Kh,1);

          gI((1:shiftInd)+oldInd,1) = [SegI_K+(nRow)*(iter-1);...
                                   SegI_Kh+(nRow)*(iter-1)];

          gJ((1:shiftInd)+oldInd,1) = [SegJ_K+(cfg.nVarSeg*2)*(iter-1);...
                                   cfg.tVar*ones((size(SegI_Kh,1)),1)];

          gV((1:shiftInd)+oldInd,1) = [SegV_K;...
                                   SegV_Kh];
%         gI( ( (1:(n1+n2)*(m1))+ ((n1+n2)*(m1+1))*(i-1)),1) = SegI+(n1+n2)*(i-1);
%         gJ( ((1:(n1+n2)*(m1))+ ((n1+n2)*(m1+1))*(i-1)),1) =  SegJ+(cfg.nVarSeg)*(i-1);
%         gV( ((1:(n1+n2)*(m1))+ ((n1+n2)*(m1+1))*(i-1)),1) =  reshape(gSegKD, [(n1+n2)*(m1),1]);
% 
%         gI( ( (1:(n1+n2))+(n1+n2)*(m1)+ ((n1+n2)*(m1+1))*(i-1)),1) = (1:(n1+n2))'+(n1+n2)*(i-1);
%         gJ( ( (1:(n1+n2))+(n1+n2)*(m1)+ ((n1+n2)*(m1+1))*(i-1)),1) = cfg.tVar*ones((n1+n2),1) ;
%         gV( ( (1:(n1+n2))+(n1+n2)*(m1)+ ((n1+n2)*(m1+1))*(i-1)),1) =  reshape(ghKD, [(n1+n2),1]);
          oldInd = oldInd+shiftInd;
          iter = iter+1;
    end
    
    g = sparse(gI(1:oldInd,1),gJ(1:oldInd,1),gV(1:oldInd,1),nRow*(cfg.nSeg-1)/2,cfg.tVar);

end

