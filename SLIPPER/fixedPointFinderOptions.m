function optParms = fixedPointFinderOptions(gVec,betaVec)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% User-defined Options
% Stored quantities: 'fval' for cost function values, 'maxlambda' for
% largest abs(eigen value)
optParms.storedQuantity='maxlambda';
optParms.useTicToc = true; % true for show calculation time

% Sampled number
optParms.sampledNumberK = 30;
optParms.sampledNumberDelta = 30;

% Range of dimensionless stiffness
optParms.kMin = 0.25;
optParms.kMax = 10;
% Range of COM velocity direction
optParms.deltaMin = -0.2;
optParms.deltaMax = -0.001;

% Range of parameters for plotting
optParms.kMinPlot = optParms.kMin;
optParms.kMaxPlot = optParms.kMax;
optParms.deltaMinPlot = optParms.deltaMin;
optParms.deltaMaxPlot = optParms.deltaMax;

%% Create parameter set for fixed-point search

if length(gVec) > 1 && length(betaVec) == 1
    optParms.searchingVarLength = length(gVec);
    optParms.g = gVec;
    optParms.beta = betaVec * ones(1, optParms.searchingVarLength);
    optParms.searchingVar = 'g';
elseif length(betaVec) > 1 && length(gVec) == 1
    optParms.searchingVarLength = length(betaVec);
    optParms.g = gVec * ones(1, optParms.searchingVarLength);
    optParms.beta = betaVec;
    optParms.searchingVar = 'beta';
elseif length(betaVec) == 1 && length(gVec) == 1
    optParms.searchingVarLength = 1;
    optParms.g = gVec;
    optParms.beta = betaVec;
    optParms.searchingVar = 'none'; % !!!!!!!!!!!!!!!!!!!!>>>>>>>>>>>>>>>>
else
    error('dimension error: only gVec or betaVec can be a vector!');
end

end

