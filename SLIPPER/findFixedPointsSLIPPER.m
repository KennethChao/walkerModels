function result = findFixedPointsSLIPPER(optParms)
%FINDFIXEDPOINTSSLIPPER Main function for finding SLIPPER fixed points
%   A unconstrained nonlinear optimization is formulated to derive the
%   fixed points of SLIPPER (SLIP with PEndulum Runner).
%
%   This function use the parfor (line 25) to speed up the calculation,
%   which requires the Parallel Computing Toolbox. To check whether the
%   toolbox is installed or not, run 'ver('distcomp')' in command line.
%   If the toolbox is not installed, please change 'parfor' to 'for' to
%   run this function.
%

%% Build buffers
[meshgridK, meshgridDelta] = createMeshGridSLIPPER(optParms);
[stablePhi, unstablePhi, stableData, unstableData, stableDataBuf, unstableDataBuf] ...
    = createDataBufferSLIPPER(optParms);
%% Opt solver option
optionsFminunc = optimset('Display', 'off', 'FinDiffType', 'central', 'MaxIter', 4e2, 'TolFun', 1e-9, 'TolX', 1e-9);

%% Find fixed points
if optParms.useTicToc
    tic
end

for k = 1:optParms.searchingVarLength
    % create buffer (used for parfor parellel computing)
    stablePhiStackBuf = nan(optParms.freeVariableNumber, optParms.sampledNumberDelta, optParms.sampledNumberK);
    unstablePhiStackBuf = nan(optParms.freeVariableNumber, optParms.sampledNumberDelta, optParms.sampledNumberK);
    
    for i = 1:optParms.sampledNumberK
        
        for j = 1:optParms.sampledNumberDelta
            parms = {};
            parms.g = optParms.g(k);
            parms.beta = optParms.beta(k);
            parms.k = meshgridK(j, i);
            
            parms.mf = optParms.mf;
            parms.rc = optParms.rc;
            parms.I = optParms.I;
            
            parms.controlMode = optParms.controlMode;
            parms.controlGain = optParms.controlGain;
            
            parms.optWeighting = optParms.optWeighting;
            
            parms.delta0 = meshgridDelta(j, i);
            
            % Initial guess of free variables
            phi0 = 0;
            phid0 = 0;
            u0 = 0;
            x0 = [phi0; phid0; u0];
            
            % Solve the unconstrained optimization
            parms.mode = 'fixedPointOpt';
            [sol, fval, exitflag, ~] = fminunc(@(x)oneStepSimulationSLIPPER(x, parms), x0, optionsFminunc);
            
            % Accept solution with positive exitflag (likely to be local
            % minimum) within cost tolerence
            if exitflag > 0 && fval < optParms.costTolerence
                
                parms.mode = 'simulationCheck';
                simResult = oneStepSimulationSLIPPER(sol, parms);
                
                % Accept solution if its duty factor is larger than 0.25
                if ~isempty(simResult.te2) && ~isempty(simResult.te) && (simResult.te / (simResult.te2 + simResult.te)) > 0.25
                    fprintf('suceed! %dth k and %dth delta\n', i, j)
                    
                    
                    % Test fixed-point stability and store sorts of results
                    parms.mode = 'perturbedSimulation';
                    PoincareMap = oneStepSimulationSLIPPER(sol, parms);
                    eigenValue = eig(PoincareMap);
                    if abs(eigenValue(1)) <= 1 && abs(eigenValue(2)) <= 1
                        stablePhiStackBuf(:, j, i) = sol;
                        stableDataBuf(j, i).fval = fval;
                        stableDataBuf(j, i).maxAbsEigenValue = max(abs(eigenValue));
                        stableDataBuf(j, i).dutyFactor = simResult.te / (simResult.te + simResult.te2);
                        stableDataBuf(j, i).runningFreqeuncy = 1 / (simResult.te + simResult.te2);
                        stableDataBuf(j, i).constantTorque = sol(3);
                        
                        springForce = -parms.k * (simResult.x(:, 1) - 1);
                        springVelocity = simResult.x(:, 2);
                        stableDataBuf(j, i).netWork = sum(springForce.*springVelocity*0.01);
                    else
                        unstablePhiStackBuf(:, j, i) = sol;
                        unstableDataBuf(j, i).fval = fval;
                        unstableDataBuf(j, i).maxAbsEigenValue = max(abs(eigenValue));
                        unstableDataBuf(j, i).dutyFactor = simResult.te / (simResult.te + simResult.te2);
                        unstableDataBuf(j, i).runningFreqeuncy = 1 / (simResult.te + simResult.te2);
                        unstableDataBuf(j, i).constantTorque = sol(3);
                        
                        springForce = -parms.k * (simResult.x(:, 1) - 1);
                        springVelocity = simResult.x(:, 2);
                        unstableDataBuf(j, i).netWork = sum(springForce.*springVelocity*0.01);
                    end
                    
                else
                    fprintf('failed! %dth k and %dth delta\n', i, j)
                end
                
            else
                fprintf('failed! %dth k and %dth delta\n', i, j)
            end
        end
    end
    % Store the result from buffer (used for parfor parellel computing)
    stablePhi(:, :, :, k) = stablePhiStackBuf;
    unstablePhi(:, :, :, k) = unstablePhiStackBuf;
    stableData(:, :, k) = stableDataBuf;
    unstableData(:, :, k) = unstableDataBuf;
end

if optParms.useTicToc
    toc
end

% Store all the results to a struct
result.meshgridK = meshgridK;
result.meshgridK = meshgridK;
result.meshgridDelta = meshgridDelta;
result.stablePhi = stablePhi;
result.unstablePhi = unstablePhi;
result.stableData = stableData;
result.unstableData = unstableData;

end
