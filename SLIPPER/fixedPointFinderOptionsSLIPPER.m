function optParms = fixedPointFinderOptionsSLIPPER(gVec,betaVec, mf, rc, I, optParms)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% User-defined Options

optParms.useTicToc = true; % true for show calculation time

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

optParms.rc = rc;
optParms.mf = mf;
optParms.I = I;

end

