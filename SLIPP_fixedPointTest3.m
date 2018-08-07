clc;
clear;
close all;
addpath('./SLIPP')
addpath('./SLIP')
warning off;

optionsFminunc = optimset('Display', 'off', 'FinDiffType', 'central', 'MaxIter', 1e4);


figure()
hold on

sampledNumber1 = 20;
sampledNumber2 = 20;

kMin = 0.25;
kMax = 20;
kMinPlot = kMin;
kMaxPlot = kMax;


deltaMin = -0.2;
deltaMax = -0.001;
deltaMinPlot = deltaMin;
deltaMaxPlot = deltaMax;


x = linspace(kMin,kMax,sampledNumber1);
y = linspace(deltaMin,deltaMax,sampledNumber2);

betaVec = 72;
% betaVec = 66:4:74;
% betaVec = 60:5:80;
% gVec = [0.025 0.05 0.1 0.21 0.46 0.66]
gVec = linspace(0.1, 0.6,6);
[X,Y] = meshgrid(x,y)
F = X;
unstableF = nan(size(X));
stableF = nan(size(X));
% for phi = 0.05:0.01:0.2
xStack = nan(2,sampledNumber2,sampledNumber1,length(betaVec));
fvalStack = nan(1,sampledNumber2,sampledNumber1);
tic

for k = 1:length(gVec)

parfor i = 1:sampledNumber1
    xStackBuf = nan(2,sampledNumber2);
    fvalStackBuf = nan(1,sampledNumber2);
    for j = 1:sampledNumber2
        
        
    parms={};
    parms.g = gVec(k);
    parms.beta = betaVec/180*pi;   
    parms.k = X(j,i);

    parms.mf = 0.5;
    parms.rc = 0.2;

%     delta0 =  -0.005;
    parms.delta0 = Y(j,i);
    phi0 =  pi/2;
    phid0 = 0;
    x0 = [phi0;phid0];

parms.mode = 'fixedPointOpt';
[x, fval, exitflag, output] = fminunc(@(x)oneStepSimulationSLIPP(x, parms), x0, optionsFminunc);

if exitflag>0 && fval<5e-4 && abs(x(1))<pi && abs(x(2))<20
    
    parms.mode = 'simulationCheck';
    ret = oneStepSimulationSLIPP(x, parms);
    if (ret.te2/ret.te)<3
    F(j,i) = fval;
    fprintf('suceed! %dth phi and %dth phid\n',i,j)
    xStackBuf(:,j) = x;
    fvalStackBuf(:,j) = fval;
    
    
    parms.mode = 'perturbedSimulation';
    ret = oneStepSimulationSLIPP(x, parms);
    e = eig(ret)
        if abs(e(1))<=1 && abs(e(2))<=1
%             stableF(j,i) = max(abs(e));
            stableF(j,i) = fval;
        else
%             unstableF(j,i) = max(abs(e));
            unstableF(j,i) = fval;
        end
    
    else
        F(j,i) = nan;    
    fprintf('failed! %dth phi and %dth phid\n',i,j)
    end
else
    F(j,i) = nan;    
    fprintf('failed! %dth phi and %dth phid\n',i,j)
end
% disp(fval)
% end
% fprintf('suceed!')
    end
  xStack(:,:,i,k)   = xStackBuf;
  fvalStack(:,:,i)   = fvalStackBuf;
end
    surf(X,Y,stableF)
    axis([kMinPlot kMaxPlot deltaMinPlot deltaMaxPlot])
end
toc