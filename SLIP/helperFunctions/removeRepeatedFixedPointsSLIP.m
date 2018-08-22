function trauncatedFixedPoints = removeRepeatedFixedPointsSLIP(fixedPoints)
%REMOVEREPEATEDFIXEDPOINTSSLIP Remove repeated fixed points from data
%   
    % Initialize 2D array "trauncatedFixedPoints"
    samplingNumbDelta = size(fixedPoints, 1);
    samplingNumbK = size(fixedPoints, 2);
    trauncatedFixedPoints = nan(samplingNumbDelta, samplingNumbK);
    
    % Search along columns of "fixedPoints", find unique numbers within in
    % the tolerence, and then store in "trauncatedFixedPoints"
    for i = 1:samplingNumbK
        buf = fixedPoints(:, i,1);
        uniqeBuf = uniquetol(buf(~isnan(buf)), 1e-3);
        if ~isempty(uniqeBuf)
            trauncatedFixedPoints(1:length(uniqeBuf), i) = uniqeBuf; %
        end
    end
    
    % Remove extra rows from "trauncatedFixedPoints"
    for i = 1:samplingNumbDelta
        if sum(isnan(trauncatedFixedPoints(i,:)))==samplingNumbK
           truncateIndex = i;
           break;
        end
    end
    trauncatedFixedPoints(truncateIndex: end, : ) = [];
end

