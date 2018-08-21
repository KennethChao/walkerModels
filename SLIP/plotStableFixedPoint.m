% close all;
clear;
clc;


%%
% h=   figure();
% set(h,'DefaultTextFontName','Liberation Serif','DefaultTextFontSize',18,...
% 'DefaultAxesFontName','Liberation Serif','DefaultAxesFontSize',18,...
% 'DefaultLineLineWidth',1,'DefaultLineMarkerSize',7.75)

% if strcmp(    storedQuantity,'delta')
%     axis([kMinFig, kMaxFig, deltaMinFig, deltaMaxFig])
% elseif( strcmp(    storedQuantity,'lambda'))
%     axis([kMinFig, kMaxFig, -1.5, 2.5])
% end
% hold on
addpath('./result/');
data = load('fixedPointData_Varing_beta_082018_2229.mat');
optParms = data.optParms
result = data.result;



for k = 1:optParms.searchingVarLength 
    for i = 1: size(result.stableSolution,1)
        yData = result.stableSolution(i,:,k);
        colorData = result.stableData(i,:,k);
        plot(optParms.kVec, yData,'k--','linewidth',0.001,'Color',[0.8,0.8,0.8]);
        scatter(optParms.kVec, yData, 12,colorData,'filled');
        hold on
    end
    
end

% legend('unstable fixed points','stable fixed points')

xlabel('$\tilde{k}$','Interpreter','latex')
ylabel('$\delta*$','Interpreter','latex')
%     ylabel('$\lambda_2$','Interpreter','latex')
axis([5, 20, 0, 0.3])
daspect([30, 1, 1])
grid on