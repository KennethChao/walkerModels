function [stablePhi,unstablePhi,stableData,unstableData] = createDataBuffer(optParms)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
stablePhi = nan(2, optParms.sampledNumberDelta, optParms.sampledNumberK, optParms.searchingVarLength);
unstablePhi = nan(2, optParms.sampledNumberDelta, optParms.sampledNumberK, optParms.searchingVarLength);
stableData = nan(optParms.sampledNumberDelta, optParms.sampledNumberK, optParms.searchingVarLength);
unstableData = nan(optParms.sampledNumberDelta, optParms.sampledNumberK, optParms.searchingVarLength);
end

