addpath('./result')
close all;
clc;
clear;

data = load('fixedPointData_Varing_none_081618_1152.mat');
optParms = data.optParms;
result = data.result;

% result.stableData = squeeze(reshape(result.stableData,1,1,[]));
    


% stableIndices = find(~isnan(result.stableData));
k=1;
    resultBuf = nan(size(result.meshgridK));
    for i = 1:optParms.sampledNumberK
        for j = 1:optParms.sampledNumberDelta
            resultBuf(j,i) = result.stableData(j,i, k).maxAbsEigenValue;
%               resultBuf(j,i) = result.stableData(j,i, k).netWork;
%               resultBuf(j,i) = result.stableData(j,i, k).dutyFactor;
%               resultBuf(j,i) = result.stableData(j,i, k).fval;
        end
    end
    
    result.meshgridK = squeeze(reshape(result.meshgridK,1,1,[]));
result.meshgridDelta = squeeze(reshape(result.meshgridDelta,1,1,[]));
result.stablePhiReshape = [squeeze(reshape(result.stablePhi(1,:,:),1,1,[])),...
                           squeeze(reshape(result.stablePhi(2,:,:),1,1,[]))];

    resultBuf = squeeze(reshape(resultBuf,1,1,[]));
stableIndices = find(~isnan(resultBuf));
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

ret = oneStepSimulationSLIPPER(x, parms);

animeSLIPPER();
end