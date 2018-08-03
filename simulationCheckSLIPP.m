addpath('./SLIPP')
close all;
x = xStack(:,1,6)
% x = [-0.1,pi/2-0.1,0]
%     parms.g = 0.12;
%     parms.beta = 72/180*pi;   
%     parms.k = 16;
    parms.g = 0.25;
    parms.beta = 72/180*pi;   
    parms.k = 4.276;
    
    parms.mf = 0.2;
    parms.rc = 0.2;
    
    parms.delta0 = -0.1;
    
%     x = [1.4788,0.0206]
parms.mode = 'simulationCheck';
% parms.mode = 'fixedPointOpt';

ret = oneStepSimulationSLIPP(x, parms)