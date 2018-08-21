addpath('./result/SLIPPER_PControl')
close all
h = figure();
hold on

colorData = 'maxAbsEigenValue';
scaleByTotalWeight = true;

% plotStableRegion('fixedPointData_Varing_none_081918_1326', colorData, scaleByTotalWeight);
plotStableRegion('fixedPointData_Varing_none_081918_1321', colorData, scaleByTotalWeight);
% plotStableRegion('fixedPointData_Varing_none_081918_1313', colorData, scaleByTotalWeight);
% plotStableRegion('fixedPointData_Varing_none_081918_1306', colorData, scaleByTotalWeight);
% plotStableRegion('fixedPointData_Varing_none_081918_1332', colorData, scaleByTotalWeight);

axis([5, 20, 0, 0.3])
daspect([20, 1, 1])
grid on


set(gcf, 'DefaultTextFontName', 'Liberation Serif', 'DefaultTextFontSize', 18, ...
    'DefaultAxesFontName', 'Liberation Serif', 'DefaultAxesFontSize', 18, ...
    'DefaultLineLineWidth', 1, 'DefaultLineMarkerSize', 7.75)


% figure labels
xlabel('$\tilde{k}$', 'Interpreter', 'latex')
ylabel('$\delta*$', 'Interpreter', 'latex')

% colorbar setup
cbh = colorbar();
titleString = 'max(abs($\lambda$))';
ylabel(cbh, titleString, 'Interpreter', 'latex')

%%

h = figure();
hold on

colorData = 'maxAbsEigenValue';
scaleByTotalWeight = true;

plotStableRegion('fixedPointData_Varing_none_081918_1321', colorData, scaleByTotalWeight);
plotStableRegion('fixedPointData_Varing_none_082018_2202', colorData, scaleByTotalWeight);
plotStableRegion('fixedPointData_Varing_none_082018_2212', colorData, scaleByTotalWeight);
plotStableRegion('fixedPointData_Varing_none_082018_2218', colorData, scaleByTotalWeight);
plotStableRegion('fixedPointData_Varing_none_082018_2222', colorData, scaleByTotalWeight);

axis([5, 20, 0, 0.3])
daspect([30, 1, 1])
grid on


set(gcf, 'DefaultTextFontName', 'Liberation Serif', 'DefaultTextFontSize', 18, ...
    'DefaultAxesFontName', 'Liberation Serif', 'DefaultAxesFontSize', 18, ...
    'DefaultLineLineWidth', 1, 'DefaultLineMarkerSize', 7.75)


% figure labels
xlabel('$\tilde{k}$', 'Interpreter', 'latex')
ylabel('$\delta*$', 'Interpreter', 'latex')

% colorbar setup
cbh = colorbar();
titleString = 'max(abs($\lambda$))';
ylabel(cbh, titleString, 'Interpreter', 'latex')
