function result = findFixedPointsSLIPPER(optParms)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

%% Build buffers
[meshgridK, meshgridDelta] = createMeshGrid(optParms);
[stablePhi, unstablePhi, stableData, unstableData] = createDataBuffer(optParms);

%% Opt solver option
optionsFminunc = optimset('Display', 'off', 'FinDiffType', 'central', 'MaxIter', 1e4);

%% Find fixed points
if optParms.useTicToc
    tic
end

for k = 1:optParms.searchingVarLength
    
    stablePhiStackBuf = nan(2, optParms.sampledNumberDelta, optParms.sampledNumberK);
    unstablePhiStackBuf = nan(2, optParms.sampledNumberDelta, optParms.sampledNumberK);
    unstableDataBuf = nan(size(meshgridK));
    stableDataBuf = nan(size(meshgridK));
    
    for i = 1:optParms.sampledNumberK
        
        parfor j = 1:optParms.sampledNumberDelta
%             for j = 1:optParms.sampledNumberDelta
            parms = {};
            parms.g = optParms.g(k);
            parms.beta = optParms.beta(k);
            parms.k = meshgridK(j, i);
            
            parms.mf = optParms.mf;
            parms.rc = optParms.rc;
            
            parms.delta0 = meshgridDelta(j, i);
            phi0 = 0;
            phid0 = 0;
            x0 = [phi0; phid0];
            
            parms.mode = 'fixedPointOpt';
            [sol, fval, exitflag, ~] = fminunc(@(x)oneStepSimulationSLIPP(x, parms), x0, optionsFminunc);
            
            if exitflag > 0 && fval < 5e-4% && abs(sol(1)) < pi && abs(sol(2)) < 20
                
                parms.mode = 'simulationCheck';
                ret = oneStepSimulationSLIPP(sol, parms);
                
                if ~isempty(ret.te2) && ~isempty(ret.te) && (ret.te2 / ret.te) < 3
                    fprintf('suceed! %dth phi and %dth phid\n', i, j)
                  
                    
                    parms.mode = 'perturbedSimulation';
                    ret = oneStepSimulationSLIPP(sol, parms);
                    eigenValue = eig(ret);
                    if abs(eigenValue(1)) <= 1 && abs(eigenValue(2)) <= 1
                        stablePhiStackBuf(:, j, i) = sol;
                        if strcmp(optParms.storedQuantity, 'fval')
                            stableDataBuf(j, i) = fval;
                        elseif strcmp(optParms.storedQuantity, 'maxlambda')
                            stableDataBuf(j, i) = max(abs(eigenValue));
                        end
                    else
                        unstablePhiStackBuf(:, j, i) = sol;
                        if strcmp(optParms.storedQuantity, 'fval')
                            unstableDataBuf(j, i) = fval;
                        elseif strcmp(optParms.storedQuantity, 'maxlambda')
                            unstableDataBuf(j, i) = max(abs(eigenValue));
                        end
                    end
                    
                else
                    fprintf('failed! %dth phi and %dth phid\n', i, j)
                end
                
            else
                fprintf('failed! %dth phi and %dth phid\n', i, j)
            end
        end
    end
    stablePhi(:, :, :, k) = stablePhiStackBuf;
    unstablePhi(:, :, :, k) = unstablePhiStackBuf;
    stableData(:, :, k) = stableDataBuf;
    unstableData(:, :, k) = unstableDataBuf;
end

if optParms.useTicToc
    toc
end

result.meshgridK = meshgridK;
result.meshgridK = meshgridK;
result.meshgridDelta = meshgridDelta;
result.stablePhi = stablePhi;
result.unstablePhi = unstablePhi;
result.stableData = stableData;
result.unstableData = unstableData;
end
