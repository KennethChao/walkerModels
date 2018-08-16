function result = findFixedPointsSLIPPER(optParms)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

%% Build buffers
[meshgridK, meshgridDelta] = createMeshGrid(optParms);
[stablePhi, unstablePhi, stableData, unstableData,stableDataBuf,unstableDataBuf] = createDataBuffer(optParms);

%% Opt solver option
optionsFminunc = optimset('Display', 'off', 'FinDiffType', 'central', 'MaxIter', 4e2,  'TolFun', 1e-12, 'TolX', 1e-12);

%% Find fixed points
if optParms.useTicToc
    tic
end

for k = 1:optParms.searchingVarLength
    
    stablePhiStackBuf = nan(2, optParms.sampledNumberDelta, optParms.sampledNumberK);
    unstablePhiStackBuf = nan(2, optParms.sampledNumberDelta, optParms.sampledNumberK);
    
    for i = 1:optParms.sampledNumberK
        
        parfor j = 1:optParms.sampledNumberDelta
            parms = {};
            parms.g = optParms.g(k);
            parms.beta = optParms.beta(k);
            parms.k = meshgridK(j, i);
            
            parms.mf = optParms.mf;
            parms.rc = optParms.rc;
            parms.I = optParms.I;
            
            parms.optWeighting = optParms.optWeighting;   
            
            parms.delta0 = meshgridDelta(j, i);
            phi0 = 0;
            phid0 = 0;
            x0 = [phi0; phid0];
            
            parms.mode = 'fixedPointOpt';
            [sol, fval, exitflag, ~] = fminunc(@(x)oneStepSimulationSLIPPER(x, parms), x0, optionsFminunc);
            
            if exitflag > 0 && fval < optParms.costTolerence
                exitflag
                parms.mode = 'simulationCheck';
                simResult = oneStepSimulationSLIPPER(sol, parms);
                
                if ~isempty(simResult.te2) && ~isempty(simResult.te) && (simResult.te/(simResult.te2 + simResult.te)) > 0.25
                    fprintf('suceed! %dth phi and %dth phid\n', i, j)
                  
                    
                    parms.mode = 'perturbedSimulation';
                    map = oneStepSimulationSLIPPER(sol, parms);
                    eigenValue = eig(map);
                    if abs(eigenValue(1)) <= 1 && abs(eigenValue(2)) <= 1
                        stablePhiStackBuf(:, j, i) = sol;
                        stableDataBuf(j, i).fval = fval;
                        stableDataBuf(j, i).maxAbsEigenValue = max(abs(eigenValue));
                        stableDataBuf(j, i).dutyFactor = simResult.te/(simResult.te+simResult.te2);
                        
                        springForce = -parms.k*(simResult.x(:,1)-1);
                        springVelocity = simResult.x(:,2);
                        stableDataBuf(j, i).netWork = sum(springForce.*springVelocity*0.01);
                    else
                        unstablePhiStackBuf(:, j, i) = sol;
                        unstableDataBuf(j, i).fval = fval;
                        unstableDataBuf(j, i).maxAbsEigenValue = max(abs(eigenValue));
                        unstableDataBuf(j, i).dutyFactor = simResult.te/(simResult.te+simResult.te2);
                        
                        springForce = -parms.k*(simResult.x(:,1)-1);
                        springVelocity = simResult.x(:,2);
                        unstableDataBuf(j, i).netWork = sum(springForce.*springVelocity*0.01);                        
                    end
                    
                else
%                     fprintf('failed! %dth phi and %dth phid\n', i, j)
                end
                
            else
%                 fprintf('failed! %dth phi and %dth phid\n', i, j)
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
