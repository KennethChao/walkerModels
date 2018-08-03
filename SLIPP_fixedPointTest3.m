clc;
clear;
close all;
addpath('./SLIPP')
addpath('./SLIP')
warning off;

optionsFminunc = optimset('Display', 'off', 'FinDiffType', 'central', 'MaxIter', 1e4);

% delta0 = 0.1;
% phi0 = pi/2+0.2;
% phid0 = -0.1;
figure()
hold on

sampledNumber1 = 30;
sampledNumber2 = 20;

x = linspace(1,20,sampledNumber1);
y = linspace(-0.3,-0.001,sampledNumber2);

betaVec = 66:4:74;

[X,Y] = meshgrid(x,y)
F = X;
% for phi = 0.05:0.01:0.2
xStack = nan(2,sampledNumber2,sampledNumber1,length(betaVec));
fvalStack = nan(1,sampledNumber2,sampledNumber1);
tic

for k = 1:length(betaVec)

parfor i = 1:sampledNumber1
    xStackBuf = nan(2,sampledNumber2);
    fvalStackBuf = nan(1,sampledNumber2);
    for j = 1:sampledNumber2
        
        
    parms={};
    parms.g = 0.2;
    parms.beta = betaVec(k)/180*pi;   
    parms.k = X(j,i);

    parms.mf = 0.2;
    parms.rc = 0.2;

    delta0 =  -0.005;
    parms.delta0 = Y(j,i);
    phi0 =  pi/2;
    phid0 = 0;
    x0 = [phi0;phid0];

parms.mode = 'fixedPointOpt';
[x, fval, exitflag, output] = fminunc(@(x)oneStepSimulationSLIPP(x, parms), x0, optionsFminunc);

if exitflag>0 && fval<5e-3 && abs(x(1))<pi && abs(x(2))<20
    
    parms.mode = 'simulationCheck';
    ret = oneStepSimulationSLIPP(x, parms);
    if (ret.te2/ret.te)<2
    
    
    F(j,i) = fval;
    fprintf('suceed! %dth phi and %dth phid\n',i,j)
    xStackBuf(:,j) = x;
    fvalStackBuf(:,j) = fval;
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
    surf(X,Y,F)
    
end
toc