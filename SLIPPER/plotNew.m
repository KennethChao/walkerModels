%%
addpath('./result')
cd ./images
clc;
close all;
clear;

data = load('fixedPointData_Varing_none_081318_1859.mat');
optParms = data.optParms;
result = data.result;

% Figure setup
h = figure();

set(h, 'DefaultTextFontName', 'Liberation Serif', 'DefaultTextFontSize', 18, ...
    'DefaultAxesFontName', 'Liberation Serif', 'DefaultAxesFontSize', 18, ...
    'DefaultLineLineWidth', 1, 'DefaultLineMarkerSize', 7.75)

% parameters for animated gif
DT = 1;

for k = 1:optParms.searchingVarLength
    C = contourf(result.meshgridK, result.meshgridDelta, result.stableData(:, :, k), 'LineStyle', ':');
%     C = surf(result.meshgridK, result.meshgridDelta, result.stableData(:, :, k), 'LineStyle', ':');
       
    if strcmp(optParms.searchingVar, 'g')
        msg = sprintf('g = %0.2f', optParms.g(k));
    elseif strcmp(optParms.searchingVar, 'beta')
        msg = sprintf('\\beta = %0.2f^o', optParms.beta(k)/pi*180);
    else
        msg = sprintf('g = %0.2f,  \\beta = %0.2f^o', optParms.g(k), optParms.beta(k)/pi*180);
    end

    % figure labels
    xlabel('$\tilde{k}$', 'Interpreter', 'latex')
    ylabel('$\delta*$', 'Interpreter', 'latex')

    % colorbar setup
    cbh = colorbar();
    titleString = 'max(abs($\lambda$))';
    ylabel(cbh, titleString, 'Interpreter', 'latex')    
    
    msg = sprintf('g = %0.2f,  \\beta = %0.2f^o', optParms.g(k), optParms.beta(k)/pi*180);
    textHandle = text(optParms.kMaxPlot*0.75, optParms.deltaMaxPlot*0.9, msg);
    msg2 = sprintf('$$\\tilde m_f = %0.2f,  \\tilde r_c = %0.2f$$',optParms.mf, optParms.rc);
    textHandle = text(optParms.kMaxPlot*0.75, optParms.deltaMaxPlot*0.75, msg2, 'Interpreter', 'latex');
    
    caxis([0.1, 1]) %round(min(min(min(result.stableData))), 1)
    axis([optParms.kMinPlot,optParms.kMaxPlot, optParms.deltaMinPlot, optParms.deltaMaxPlot])
    axis([10,20, optParms.deltaMinPlot, optParms.deltaMaxPlot])
    axis([10,15, optParms.deltaMinPlot, 0.1])
    
    fileName = sprintf('fixedPointData_Varing_%s_0%d.tif', optParms.searchingVar, k);
    saveas(h,fileName);
    pause(1.0)
end
cd ../