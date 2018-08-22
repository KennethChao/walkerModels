function reshapedFixedPoints  = reshapeFixedPointsSLIP(stableFixedPoints,unstableFixedPoints,kVec,deltaMax,startingPointOption)
%RESHAPEFIXEDPOINTSSLIP Sorting the data based on the speicified starting
% point option
%  
    


    blockLength = (size(stableFixedPoints,1)+size(unstableFixedPoints,1));
    
    bufferStableFixedPoints = reshape(stableFixedPoints',1,[]);
    bufferUnstableFixedPoints = reshape(unstableFixedPoints',1,[]);
    
    bufferFixedPoints = [bufferStableFixedPoints,bufferUnstableFixedPoints];
    reshapeFactor = max(kVec) /deltaMax;    
    bufferKVec = zeros(1,blockLength*length(kVec));
    for i = 1:blockLength
        bufferKVec(( 1:length(kVec) )+ length(kVec)*(i-1)  ) = kVec;
    end
    
    nanIndex = find(isnan(bufferFixedPoints));
    bufferFixedPoints(nanIndex) = [];
    bufferKVec(nanIndex) = [];
    

    
    bufferData = [bufferFixedPoints*reshapeFactor;bufferKVec];
    bufferDataNorm = zeros(1,size(bufferData,2));
    
    reshapedFixedPoints = nan(size(bufferData));
    
    for i = 1:size(bufferData,2)
        bufferDataNorm(i) = norm(bufferData(:,i));
    end
    if strcmp(startingPointOption,'minNorm')
    [~,normIndex] = min(bufferDataNorm);
    elseif strcmp(startingPointOption,'maxNorm')
    [~,normIndex] = max(bufferDataNorm);        
    elseif strcmp(startingPointOption,'minDelta')
    [~,normIndex] = min(bufferData(1,:));        
    end
    
    reshapedFixedPoints(:,1) = bufferData(:,normIndex);
    bufferData(:,normIndex) = [];
    
    for i = 1:size(bufferData,2)
        bufferDataNorm = zeros(1,size(bufferData,2));
        for j = 1:size(bufferData,2)
            bufferDataNorm(j) = norm(bufferData(:,j)-reshapedFixedPoints(:,i));
        end
        [~,normIndex] = min(bufferDataNorm);
    
        reshapedFixedPoints(:,i+1) = bufferData(:,normIndex);
        bufferData(:,normIndex) = [];   
    end

    reshapedFixedPoints(1,:) = reshapedFixedPoints(1,:)/reshapeFactor;


end

