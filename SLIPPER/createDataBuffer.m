function [stablePhi,unstablePhi,stableData,unstableData,stableDataBuf,unstableDataBuf] = createDataBuffer(optParms)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
stablePhi = nan(2, optParms.sampledNumberDelta, optParms.sampledNumberK, optParms.searchingVarLength);
unstablePhi = nan(2, optParms.sampledNumberDelta, optParms.sampledNumberK, optParms.searchingVarLength);
% stableData = nan(optParms.sampledNumberDelta, optParms.sampledNumberK, optParms.searchingVarLength);
% unstableData = nan(optParms.sampledNumberDelta, optParms.sampledNumberK, optParms.searchingVarLength);


for i = 1:optParms.sampledNumberK
    for j = 1:optParms.sampledNumberDelta
        unstableDataBuf(j,i).fval = nan;
        unstableDataBuf(j, i).maxAbsEigenValue = nan;
        unstableDataBuf(j, i).dutyFactor = nan;
        unstableDataBuf(j, i).netWork = nan;
        stableDataBuf(j, i).fval = nan;
        stableDataBuf(j, i).maxAbsEigenValue = nan;
        stableDataBuf(j, i).dutyFactor = nan;
        stableDataBuf(j, i).netWork = nan;
    end
end

for i = 1:optParms.sampledNumberK
    for j = 1:optParms.sampledNumberDelta
        for k = 1:optParms.searchingVarLength
            unstableData(j,i,k).fval = nan;
            unstableData(j, i,k).maxAbsEigenValue = nan;
            unstableData(j, i,k).dutyFactor = nan;
            unstableData(j, i,k).netWork = nan;
            stableData(j, i,k).fval = nan;
            stableData(j, i,k).maxAbsEigenValue = nan;
            stableData(j, i,k).dutyFactor = nan;
            stableData(j, i,k).netWork = nan;
        end
    end
end


end

