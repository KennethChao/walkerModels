addpath('./result')
close all;
clc;
clear;

data = load('fixedPointData_Varing_none_081018_1151.mat');
optParms = data.optParms;
result = data.result;

result.meshgridK = squeeze(reshape(result.meshgridK,1,1,[]));
result.meshgridDelta = squeeze(reshape(result.meshgridDelta,1,1,[]));
result.stablePhiReshape = [squeeze(reshape(result.stablePhi(1,:,:),1,1,[])),...
                           squeeze(reshape(result.stablePhi(2,:,:),1,1,[]))];
result.stableData = squeeze(reshape(result.stableData,1,1,[]));
    


stableIndices = find(~isnan(result.stableData));
%     x = [1.4788,0.0206]


%%

% testIndex =12;
% stableIndex = stableIndices(testIndex);

for i = 1:length(stableIndices)
stableIndex = stableIndices(i);
% fprintf('max abs(eigen value):%.3f\n',result.stableData(stableIndex));
parms = optParms;
parms.mode = 'simulationCheck';
parms.k = result.meshgridK(stableIndex);
parms.delta0 = result.meshgridDelta(stableIndex);

x = result.stablePhiReshape(stableIndex,:);
% fprintf('phi = %.3f, dot phi = %.3f\n',x(1),x(2));

ret = oneStepSimulationSLIPP(x, parms);

animeSLIPP();
end