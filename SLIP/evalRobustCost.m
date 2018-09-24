% function plotAllFixedPointsSLIP(filename,dataType,startingPointOption)
% PLOTALLFIXEDPOINTSSLIP helper function to plot all fixed-point set of SLIP
% (for comparison between SLIP and SLIPPER)
%
close all
% Load data
addpath('./helperFunctions')
addpath('./result')
data = load('fixedPointData_Varing_none_091418_1109.mat');
dataType = 'delta';

optParms = data.optParms;
disp(optParms);
result = data.result;

for k = 1:optParms.searchingVarLength
    
    % Remove repeated fixed point solutions
    if strcmp(dataType,'eigenValue')
        trauncatedUnstableSolution = removeRepeatedFixedPointsSLIP(result.unstableData(:, :,k));
        trauncatedStableSolution = removeRepeatedFixedPointsSLIP(result.stableData(:, :,k));    
    elseif strcmp(dataType,'delta')
        trauncatedUnstableSolution = removeRepeatedFixedPointsSLIP(result.unstableSolution(:, :,k));
        trauncatedUnstableData = removeRepeatedFixedPointsSLIP(result.unstableData(:, :,k));    
        trauncatedStableSolution = removeRepeatedFixedPointsSLIP(result.stableSolution(:, :,k));        
        trauncatedStableData = removeRepeatedFixedPointsSLIP(result.stableData(:, :,k));    
    end

%     hold on
    index = isnan(trauncatedStableSolution);
    trauncatedStableSolution(index)=[];
    optParms.kVec(index) = [];
    trauncatedStableData(index) = [];
end



iterationNumber = length(trauncatedStableSolution);

robustCost = nan(size(trauncatedStableSolution));
data = load('perturbationABuffer_0914.mat','disturbanceBuffer');
robustCostBuffer = [];
disturbanceBuffer = data.disturbanceBuffer;
surviveStepsBuffer = [];
for k = 1:10

for i = 1:iterationNumber
    delta0 = trauncatedStableSolution(i);
    optParms.k = optParms.kVec(i);
    optParms.mode = 'dataCollection';
    
    result = oneStepSimulationSLIP(delta0, optParms);
    
    totalTime = result.te + result.te2
    
    tInterval = 0.15;
    
    timeStack = linspace(totalTime-tInterval, totalTime+tInterval,20) - result.te;
    
    cost = 0;
    for j = 1:length(timeStack)
        optParms.mode = 'fixedPointOpt';
        t2 = timeStack(j);
        result = oneStepSimulationSLIP_SpecifiedTime(delta0,t2, optParms);
        cost = cost +result;
    end
    robustCost(i)=cost;
end
robustCostBuffer = [robustCostBuffer;robustCost];

% disturbance = 0.08*rand(1,200)-0.04;
disturbance = disturbanceBuffer(k,:)
% disturbanceBuffer = [disturbanceBuffer;disturbance];

surviveSteps = nan(size(trauncatedStableSolution));
for i = 1:iterationNumber
    delta0 = trauncatedStableSolution(i);
    optParms.k = optParms.kVec(i);
    optParms.mode = 'robustnessTest';
    
    result = oneStepPerturbedSimulationSLIP(delta0,disturbance, optParms);
    surviveSteps(i) = result;
%     totalTime = result.te + result.te2
    
%     tInterval = 0.05*totalTime;
    
%     timeStack = linspace(totalTime-tInterval, totalTime+tInterval,10) - result.te;
    
%     cost = 0;
%     for j = 1:length(timeStack)
%         optParms.mode = 'fixedPointOpt';
%         t2 = timeStack(j);
%         result = oneStepSimulationSLIP_SpecifiedTime(delta0,t2, optParms);
%         cost = cost +result;
%     end
%     robustCost(i)=cost/tInterval;
    surviveStepsBuffer = [surviveStepsBuffer;surviveSteps];
end

end
figure()
plot(mean(robustCostBuffer))
figure()
plot(nanmean(surviveStepsBuffer))