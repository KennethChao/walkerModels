function c = constDym(x,dx,ddx,h,parms)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


dim = size(x,1);
phaseNum = length(parms.phase);

c = zeros(1,parms.totalKnotNumber*dim);
shiftIndex = 0;

for i = 1:phaseNum
    for j =  1: (parms.phase(i).KnotNumber)
    
    cSegment = ddx(:,j) - parms.phase(i).dymFunc(x(:,j),dx(:,j), parms);
        
    c(shiftIndex + (1:dim)) = cSegment;
    shiftIndex = shiftIndex + dim;
    end

end

