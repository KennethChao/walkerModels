addpath('./result')
addpath('./result/SLIPPER_PControl')
close all;
clc;
clear;

data = load('fixedPointData_Varing_none_082018_1403.mat');
optParms = data.optParms
result = data.result;


for k = 1:1 %1:optParms.searchingVarLength
    result = data.result;
    resultBuf = nan(size(result.meshgridK));
    for i = 1:optParms.sampledNumberK
        for j = 1:optParms.sampledNumberDelta
            resultBuf(j, i) = result.stableData(j, i, k).maxAbsEigenValue;
        end
    end
    stableIndices = find(~isnan(resultBuf));
    
    
    for i = 1:optParms.sampledNumberK
        for j = 1:optParms.sampledNumberDelta
            resultBuf(j, i) = result.stableData(j, i, k).maxAbsEigenValue;
            %             resultBuf(j,i) = result.stableData(j,i, k).netWork;
            %               resultBuf(j,i) = result.stableData(j,i, k).runningFreqeuncy;
            %               resultBuf(j,i) = result.stableData(j,i, k).dutyFactor;
            %               resultBuf(j,i) = result.stableData(j,i, k).constantTorque;
            %               resultBuf(j,i) = result.stableData(j,i, k).fval;
            %               resultBuf(j,i) = result.meshgridK(j,i);
            %               resultBuf(j,i) = result.meshgridDelta(j,i);
        end
    end
    resultBuf = squeeze(reshape(resultBuf, 1, 1, []));
    
    result.meshgridK = squeeze(reshape(result.meshgridK, 1, 1, []));
    result.meshgridDelta = squeeze(reshape(result.meshgridDelta, 1, 1, []));
    result.stablePhiReshape = [squeeze(reshape(result.stablePhi(1, :, :), 1, 1, [])), ...
        squeeze(reshape(result.stablePhi(2, :, :), 1, 1, [])), ...
        squeeze(reshape(result.stablePhi(3, :, :), 1, 1, []))];
    
    
    plotParms.stableFixedPointNumber = length(stableIndices);
    plotParms.colorMapMax = max(resultBuf);
    plotParms.colorMapMin = min(resultBuf);
    
    %%
    for i = 1:length(stableIndices)
        stableIndex = stableIndices(i);
        parms = optParms;
        parms.mode = 'simulationCheck';
        parms.k = result.meshgridK(stableIndex);
        parms.delta0 = result.meshgridDelta(stableIndex);
        
        parms.beta = optParms.beta(k);
        parms.g = optParms.g(k);
        
        parms.controlMode = optParms.controlMode;
        parms.controlGain = optParms.controlGain;
        
        parms.optWeighting = optParms.optWeighting;
        
        x = result.stablePhiReshape(stableIndex, :);
        ret = oneStepSimulationSLIPPER(x, parms);
        
        %         'phasePotrait_zf', 'phasePotrait_phiVec&zf', 'phasePotrait_zVec&phi'
                phasePortraitSLIPPER('phasePotrait_phi',ret,resultBuf(stableIndex),parms,plotParms);
%         phasePortraitSLIPPER('phasePotrait_phiVec&zf', ret, resultBuf(stableIndex), parms, plotParms);
%                 phasePortraitSLIPPER('phasePotrait_zVec&phi',ret,resultBuf(stableIndex),parms,plotParms);
        
    end
    
end

cbh = colorbar;