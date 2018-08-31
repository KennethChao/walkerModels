function c = constKineHSM(x,dx,ddx,h,parms)
%CONSTKINEHSM Hermite-Simpson constraints
%   Detailed explanation goes here

dim = size(x,1)*2;

for i = 1: parms.phaseNum
    if parms.phase(i).knotNumber<=2
        error('knotNumber is too small!');
    end
end

c = zeros(1,parms.totaHSMCnstNumber*dim*2);

iteration = 0;

for i = 1:parms.phaseNum
    if i ==1
        indexRange = 1: parms.phase(i).knotNumber;
        state = [x(:, indexRange);dx(:,indexRange)];
        dState = [dx(:,indexRange);ddx(:,indexRange)];            
    else        
        indexRange = parms.phase(i-1).knotNumber + (1: parms.phase(i).knotNumber);
        state = [x(:, indexRange);dx(:,indexRange)];
        dState = [dx(:,indexRange);ddx(:,indexRange)];                    
    end
    for j =  1:2: (parms.phase(i).knotNumber-2)
        cSegment1 = state(:,j+2)-state(:,j)-1/6*h(i)*(dState(:,j) + 4*dState(:,j+1)+dState(:,j+2));
        cSegment2 = state(:,j+1)-1/2*(state(:,j) + state(:,j+2))-1/8*h(i)*(dState(:,j)-dState(:,j+2));

%         if i ==1
%         indexRange = 2*dim*(iteration) + (1:^2*dim);
%         else
        indexRange = 2*dim*(iteration) + (1:2*dim);% + (parms.phase(i-1).knotNumber-1)*2*dim;            
%         end
        
        iteration = iteration + 1;
        c(1,indexRange) = [cSegment1',cSegment2'];
    end

end

