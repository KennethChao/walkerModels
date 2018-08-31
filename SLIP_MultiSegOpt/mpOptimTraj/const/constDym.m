function c = constDym(x,dx,ddx,parms)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


dim = size(x,1);
phaseNum = length(parms.phase);

c = zeros(1,parms.totalKnotNumber*dim);
shiftIndex = 0;

for i = 1:phaseNum
    
    if i ==1
        indexRange = 1: parms.phase(i).knotNumber;        
        xSeg = x(:,indexRange);
        dxSeg = dx(:,indexRange);      
        ddxSeg = ddx(:,indexRange);      
    else        
        indexRange = parms.phase(i-1).knotNumber + (1: parms.phase(i).knotNumber);   
        xSeg = x(:,indexRange);
        dxSeg = dx(:,indexRange);        
        ddxSeg = ddx(:,indexRange);      
    end       
    
    for j =  1: (parms.phase(i).knotNumber)
    
    cSegment = ddxSeg(:,j) - parms.phase(i).dymFunc(xSeg(:,j),dxSeg(:,j), parms.g, parms.k);
        
    c(shiftIndex + (1:dim)) = cSegment;
    shiftIndex = shiftIndex + dim;
    end

end

