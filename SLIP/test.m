  parms={};
        parms.g = 0.46;
        parms.beta = 72/180*pi; 
        parms.k = 16.03;
%     parms.mode = 'perturbedSimulation';
    parms.mode = 'fixedPointOpt';
result = oneStepSimulationSLIP(0.09331, parms)