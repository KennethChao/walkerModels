%% Plot animated gif
addpath('./retsult')
cd ./images
clc;
close all;
clear;


data = load('fixedPointData_Varing_none_080918_1054.mat');
optParms = data.optParms;
result = data.result;

% Figure setup
h = figure();

set(h, 'DefaultTextFontName', 'Liberation Serif', 'DefaultTextFontSize', 18, ...
    'DefaultAxesFontName', 'Liberation Serif', 'DefaultAxesFontSize', 18, ...
    'DefaultLineLineWidth', 1, 'DefaultLineMarkerSize', 7.75)
hold on

% figure labels
xlabel('$\tilde{k}$', 'Interpreter', 'latex')
ylabel('$\delta*$', 'Interpreter', 'latex')

% colorbar setup
cbh = colorbar();
titleString = 'max(abs($\lambda$))';
ylabel(cbh, titleString, 'Interpreter', 'latex')

% text for chaging parameters
if strcmp(optParms.searchingVar, 'g')
    msg = sprintf('g = %0.2f', optParms.g(1));
elseif strcmp(optParms.searchingVar, 'beta')
    msg = sprintf('beta = %0.2^of', optParms.beta(1));
end
textHandle = text(optParms.kMaxPlot/2, optParms.deltaMinPlot/3, msg);

% parameters for animated gif
DT = 1;

% fileName
if strcmp(optParms.searchingVar, 'g')
    filename = 'SLIPPER_FixedPointAnimation_VaryingG.gif'; % Specify the output file name
else
    filename = 'SLIPPER_FixedPointAnimation_VaryingBeta.gif'; % Specify the output file name
end

% Start plotting
for k = 1:optParms.searchingVarLength
    C = contourf(result.meshgridK, result.meshgridDelta, result.stableData(:, :, k), 'LineStyle', ':');
       
    if strcmp(optParms.searchingVar, 'g')
        msg = sprintf('g = %0.2f', optParms.g(k));
    elseif strcmp(optParms.searchingVar, 'beta')
        msg = sprintf('beta = %0.2f^o', optParms.beta(k)/pi*180);
    end
    
    set(textHandle, 'String', msg)
    caxis([0.8, 1]) %round(min(min(min(result.stableData))), 1)
    axis([optParms.kMinPlot,optParms.kMaxPlot, optParms.deltaMinPlot, optParms.deltaMaxPlot])
    
    frame = getframe(h);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    if k == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', DT);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', DT);
    end
    pause(1.0)
end

cd ../../
%%
addpath('./SLIPPER/ret')
cd ./SLIPPER/images
clc;
close all;
clear;

data = load('fixedPointData_Varing_beta_080818_1239.mat');
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
       
    if strcmp(optParms.searchingVar, 'g')
        msg = sprintf('g = %0.2f', optParms.g(k));
    elseif strcmp(optParms.searchingVar, 'beta')
        msg = sprintf('\\beta = %0.2f^o', optParms.beta(k)/pi*180);
    end

    % figure labels
    xlabel('$\tilde{k}$', 'Interpreter', 'latex')
    ylabel('$\delta*$', 'Interpreter', 'latex')

    % colorbar setup
    cbh = colorbar();
    titleString = 'max(abs($\lambda$))';
    ylabel(cbh, titleString, 'Interpreter', 'latex')    
    
    textHandle = text(optParms.kMaxPlot/2, optParms.deltaMinPlot/4, msg);
    caxis([0.8, 1]) %round(min(min(min(result.stableData))), 1)
    axis([optParms.kMinPlot,optParms.kMaxPlot, optParms.deltaMinPlot, optParms.deltaMaxPlot])
    
    fileName = sprintf('fixedPointData_Varing_%s_0%d.tif', optParms.searchingVar, k);
    saveas(h,fileName);
    pause(1.0)
end
cd ../../