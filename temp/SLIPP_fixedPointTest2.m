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

sampledNumber1 = 5;
sampledNumber2 = 5;

x = linspace(pi/2-0.5,pi/2,sampledNumber1);
y = linspace(-0.5,0.5,sampledNumber2);

[X,Y] = meshgrid(x,y)
F = X;
% for phi = 0.05:0.01:0.2
xStack = nan(2,sampledNumber2,sampledNumber1);
fvalStack = nan(1,sampledNumber2,sampledNumber1);
tic
parfor i = 1:sampledNumber1
    xStackBuf = nan(2,sampledNumber2);
    fvalStackBuf = nan(1,sampledNumber2);
    for j = 1:sampledNumber2
        
        
    parms={};
    parms.g = 0.46;
    parms.beta = 72/180*pi;   
    parms.k = 16;

    parms.mf = 0.2;
    parms.rc = 0.2;

    delta0 =  -0.005;
    parms.delta0 = delta0;
    phi0 =  X(j,i);
    
    phid0 = Y(j,i);
    x0 = [phi0;phid0];

parms.mode = 'fixedPointOpt';
[x, fval, exitflag, output] = fminunc(@(x)oneStepSimulationSLIPP(x, parms), x0, optionsFminunc);

if exitflag>0 && fval<1e-2 && abs(x(1))<pi && abs(x(2))<20
    F(j,i) = fval;
    fprintf('suceed! %dth phi and %dth phid\n',i,j)
    xStackBuf(:,j) = x;
    fvalStackBuf(:,j) = fval;
else
    F(j,i) = nan;    
    fprintf('failed! %dth phi and %dth phid\n',i,j)
end
% disp(fval)
% end
% fprintf('suceed!')
    end
  xStack(:,:,i)   = xStackBuf;
  fvalStack(:,:,i)   = fvalStackBuf;
end
toc
% %   x0         = [delta0 phi0 phid0];  % The starting point.
% %   options.cl = [ 0     0    -1 ];             % Lower bounds on constraints.
% %   options.cu = [ pi/2  pi   1 ];             % Upper bounds on constraints.
% % 
% %   % Set the IPOPT options.
% % %   options.ipopt.jac_c_constant        = 'yes';
% %   options.ipopt.hessian_approximation = 'limited-memory';
% %   options.ipopt.jacobian_approximation = 'finite-difference-values';
% % %   options.ipopt.mu_strategy           = 'adaptive';
% % %   options.ipopt.tol                   = 1e-7;
% % 
% %   % The callback functions.
% %   funcs.objective         = @cost;
% %   funcs.gradient          = @gcost;
% % %   funcs.constraints       = @(x)oneStepSimulationSLIPP(x, parms);
% % %   funcs.jacobian          = @(x) sparse(ones(1,3));
% % %   funcs.jacobianstructure = @(x) sparse(ones(1,3));
% %   % Run IPOPT.
% %   [x info] = ipopt(x0,funcs,options);