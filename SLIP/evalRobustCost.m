% function plotAllFixedPointsSLIP(filename,dataType,startingPointOption)
% PLOTALLFIXEDPOINTSSLIP helper function to plot all fixed-point set of SLIP
% (for comparison between SLIP and SLIPPER)
%

% Load data
data = load('fixedPointData_Varing_none_091118_1519.mat');
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





for i = 1:iterationNumber
    delta0 = trauncatedStableSolution(i);
    optParms.k = optParms.kVec(i);
    optParms.mode = 'dataCollection';
    
    result = oneStepSimulationSLIP(delta0, optParms);
    
    totalTime = result.te + result.te2
    
    tInterval = 0.2*totalTime;
    
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
plot(robustCost)

disturbance = 0.1*rand(1,200)-0.05;

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
end

