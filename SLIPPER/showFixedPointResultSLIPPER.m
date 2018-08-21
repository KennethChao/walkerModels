% This is the script to plot fixed-point set of SLIPPER
% (for comparison between SLIP and SLIPPER)
% 
% In the SLIPPER model, an inertia and a p control for phi dot are
% included. If 'scaleByTotalMass' is true, the plot will scaled so that
% the total mass will scaled to 1.

addpath('./result/SLIPPER_PControl')
addpath('./helperFunctions')
close all



%% Fixed beta, varying dimensionless g
% plot the fixed points in the order of:
% beta = 72 degree, dimensionless g = [0.12, 0.2, 0.28, 0.36, 0.44]
h1 = figure();
hold on

colorData = 'maxAbsEigenValue';
scaleByTotalMass = true;

plotStableRegionSLIPPER('fixedPointData_Varing_none_081918_1326', colorData, scaleByTotalMass);
plotStableRegionSLIPPER('fixedPointData_Varing_none_081918_1321', colorData, scaleByTotalMass);
plotStableRegionSLIPPER('fixedPointData_Varing_none_081918_1313', colorData, scaleByTotalMass);
plotStableRegionSLIPPER('fixedPointData_Varing_none_081918_1306', colorData, scaleByTotalMass);
plotStableRegionSLIPPER('fixedPointData_Varing_none_081918_1332', colorData, scaleByTotalMass);

axis([5, 20, 0, 0.3])
daspect([20, 1, 1])
grid on

% font setup
set(gcf, 'DefaultTextFontName', 'Liberation Serif', 'DefaultTextFontSize', 18, ...
    'DefaultAxesFontName', 'Liberation Serif', 'DefaultAxesFontSize', 18, ...
    'DefaultLineLineWidth', 1, 'DefaultLineMarkerSize', 7.75)


% figure labels
xlabel('$\tilde{k}$', 'Interpreter', 'latex')
ylabel('$\delta^*$', 'Interpreter', 'latex')

% colorbar setup
cbh = colorbar();
titleString = 'max(abs($\lambda$))';
ylabel(cbh, titleString, 'Interpreter', 'latex')

%% Fixed dimensionless g, varying beta
% plot the fixed points in the order of:
% dimensionless g = 0.2, beta = [66 68 70 72 74] (degree)
h2 = figure();
hold on

colorData = 'maxAbsEigenValue';
scaleByTotalMass = true;


plotStableRegionSLIPPER('fixedPointData_Varing_none_082018_2202', colorData, scaleByTotalMass);
plotStableRegionSLIPPER('fixedPointData_Varing_none_081918_1321', colorData, scaleByTotalMass);
plotStableRegionSLIPPER('fixedPointData_Varing_none_082018_2212', colorData, scaleByTotalMass);
plotStableRegionSLIPPER('fixedPointData_Varing_none_082018_2218', colorData, scaleByTotalMass);
plotStableRegionSLIPPER('fixedPointData_Varing_none_082018_2222', colorData, scaleByTotalMass);

axis([5, 20, 0, 0.3])
daspect([30, 1, 1])
grid on

% font setup
set(gcf, 'DefaultTextFontName', 'Liberation Serif', 'DefaultTextFontSize', 18, ...
    'DefaultAxesFontName', 'Liberation Serif', 'DefaultAxesFontSize', 18, ...
    'DefaultLineLineWidth', 1, 'DefaultLineMarkerSize', 7.75)

% figure labels
xlabel('$\tilde{k}$', 'Interpreter', 'latex')
ylabel('$\delta^*$', 'Interpreter', 'latex')

% colorbar setup
cbh = colorbar();
titleString = 'max(abs($\lambda$))';
ylabel(cbh, titleString, 'Interpreter', 'latex')
