function [stablePhi, unstablePhi, stableData, unstableData, stableDataBuf, unstableDataBuf] = createDataBuffer(optParms)
%CREATEDATABUFFER Create buffer for storing results of SLIPPER fixed points
%   Create double arrays for stable and unstable fixed points, and 
%   strucutre arrays of structure for related results.
%
stablePhi = nan(optParms.freeVariableNumber, optParms.sampledNumberDelta, optParms.sampledNumberK, optParms.searchingVarLength);
unstablePhi = nan(optParms.freeVariableNumber, optParms.sampledNumberDelta, optParms.sampledNumberK, optParms.searchingVarLength);

for i = 1:optParms.sampledNumberK
    for j = 1:optParms.sampledNumberDelta
        unstableDataBuf(j, i).fval = nan;
        unstableDataBuf(j, i).maxAbsEigenValue = nan;
        unstableDataBuf(j, i).dutyFactor = nan;
        unstableDataBuf(j, i).netWork = nan;
        unstableDataBuf(j, i).runningFreqeuncy = nan;
        unstableDataBuf(j, i).constantTorque = nan;
        stableDataBuf(j, i).fval = nan;
        stableDataBuf(j, i).maxAbsEigenValue = nan;
        stableDataBuf(j, i).dutyFactor = nan;
        stableDataBuf(j, i).netWork = nan;
        stableDataBuf(j, i).runningFreqeuncy = nan;
        stableDataBuf(j, i).constantTorque = nan;
    end
end

for i = 1:optParms.sampledNumberK
    for j = 1:optParms.sampledNumberDelta
        for k = 1:optParms.searchingVarLength
            unstableData(j, i, k).fval = nan;
            unstableData(j, i, k).maxAbsEigenValue = nan;
            unstableData(j, i, k).dutyFactor = nan;
            unstableData(j, i, k).netWork = nan;
            unstableData(j, i).runningFreqeuncy = nan;
            unstableData(j, i).constantTorque = nan;
            stableData(j, i, k).fval = nan;
            stableData(j, i, k).maxAbsEigenValue = nan;
            stableData(j, i, k).dutyFactor = nan;
            stableData(j, i, k).netWork = nan;
            stableData(j, i).runningFreqeuncy = nan;
            stableData(j, i).constantTorque = nan;
        end
    end
end


end
