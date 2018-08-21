close all;
clc;


%%
% h=   figure();
% set(h,'DefaultTextFontName','Liberation Serif','DefaultTextFontSize',18,...
% 'DefaultAxesFontName','Liberation Serif','DefaultAxesFontSize',18,...
% 'DefaultLineLineWidth',1,'DefaultLineMarkerSize',7.75)

if strcmp(    storedQuantity,'delta')
    axis([kMinFig, kMaxFig, deltaMinFig, deltaMaxFig])
elseif( strcmp(    storedQuantity,'lambda'))
    axis([kMinFig, kMaxFig, -1.5, 2.5])
end
hold on

for k = 1:5
    
    trauncatedUnstableSolution = removeRepeatedFixedPoints(unstableSolution(:, :,k));
    trauncatedStableSolution = removeRepeatedFixedPoints(stableSolution(:, :,k));
    
%     plot(kVec, trauncatedUnstableSolution(1, :), 'ro','MarkerSize',5,'MarkerFaceColor' , [1 0 0] );
%     plot(kVec, trauncatedStableSolution(1, :), 'bo','MarkerSize',5,'MarkerFaceColor' , [0 0 1] );
    
%     for i = 1: size(trauncatedUnstableSolution,1)
%         h1 = plot(kVec, trauncatedUnstableSolution(i, :), 'r-','MarkerSize',5,'MarkerFaceColor' , [1 0 0] );
%     end
    for i = 1: size(trauncatedStableSolution,1)
        h2 = plot(kVec, trauncatedStableSolution(i, :), 'b-','MarkerSize',1,'MarkerFaceColor' , [0 0 1] );
    end

  startingPointOption = 'minDelta';
%   startingPointOption = 'minNorm2Zero';
%   startingPointOption = 'maxNorm2Zero';
  reshapedFixedPoints  = reshapeFixedPoints(trauncatedStableSolution,trauncatedUnstableSolution,kVec,deltaMax,startingPointOption);

% ax = gca;
% ax.ColorOrderIndex = 1;

%     plot(reshapedFixedPoints(2,:),reshapedFixedPoints(1,:))  
end

% legend('unstable fixed points','stable fixed points')

xlabel('$\tilde{k}$','Interpreter','latex')
ylabel('$\delta*$','Interpreter','latex')
%     ylabel('$\lambda_2$','Interpreter','latex')
