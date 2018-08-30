function c = constKineHSM(x,dx,ddx,h,parms)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


dim = size(x,1)*2;
phaseNum = length(parms.phase);

for i = 1: phaseNum
    if parms.phase(i).knotNumber<=2
        error('knotNumber is too small!');
    end
end


c = zeros(1,(parms.totalKnotNumber-4)*dim);
shiftIndex = 0;

for i = 1:phaseNum
    for j =  1: (parms.phase(i).knotNumber-2)

    state = [x(:,shiftIndex+(1:parms.phase(i).knotNumber));dx(:,shiftIndex+(1:parms.phase(i).knotNumber))];
    dState = [dx(:,shiftIndex+(1:parms.phase(i).knotNumber));ddx(:,shiftIndex+(1:parms.phase(i).knotNumber))];
    
    cSegment1 = state(:,j+2)-state(:,j)-1/6*h(i)*(dState(:,j) + 4*dState(:,j+1)+dState(:,j+2));
    cSegment2 = state(:,j+1)-1/2*(state(:,j) + state(:,j+2))-1/8*h(i)*(dState(:,j)-dState(:,j+2));
    
    c(1,2*dim*(i-1) + (1:2*dim)) = [cSegment1',cSegment2'];
    shiftIndex = shiftIndex + parms.phase(i).knotNumber;
    end

end

