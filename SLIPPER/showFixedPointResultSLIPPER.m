% This is the script to plot stable fixed-point set of SLIPPER
% (for comparison between SLIP and SLIPPER)
%
% In the SLIPPER model, an inertia and a p control for phi dot are
% included. If 'scaleByTotalMass' is true, the plot will scaled so that
% the total mass will scaled to 1.

addpath('./result/SLIPPER_PControl')
addpath('./helperFunctions')

addpath('../SLIP/result/Comparison2SLIPPER')
addpath('../SLIP/helperFunctions')

close all

%% Delta-dimensiolness k plot: Fixed beta, varying dimensionless g
% plot the fixed points in the order of:
% beta = 72 degree, dimensionless g = [0.12, 0.2, 0.28, 0.36, 0.44]
h1 = figure();
hold on

colorData = 'maxAbsEigenValue';
scaleByTotalMass = true;

plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_082218_1529', colorData, scaleByTotalMass);
plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_082218_1604', colorData, scaleByTotalMass);
plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_082218_1617', colorData, scaleByTotalMass);
plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_082218_1644', colorData, scaleByTotalMass);
plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_082218_1649', colorData, scaleByTotalMass);

plotStableFixedPointsSLIP('fixedPointData_Varing_g_dutyFactorChecked_082118_1954');

axis([5, 20, 0, 0.3])
daspect([30, 1, 1])
grid on

% font setup
set(gcf, 'DefaultTextFontName', 'Liberation Serif', 'DefaultTextFontSize', 18, ...
    'DefaultAxesFontName', 'Liberation Serif', 'DefaultAxesFontSize', 18, ...
    'DefaultLineLineWidth', 1, 'DefaultLineMarkerSize', 7.75)

% add parameter text
text(8,0.32,'$\beta = 72^o, \tilde{g} = [0.12, 0.2, 0.28, 0.36, 0.44]$','Interpreter','latex','FontSize',12)


% figure labels
xlabel('$\tilde{k}$', 'Interpreter', 'latex')
ylabel('$\delta^*$', 'Interpreter', 'latex')

% colorbar setup
cbh = colorbar();
titleString = 'max(abs($\lambda$))';
ylabel(cbh, titleString, 'Interpreter', 'latex')

%% Phi-phid plot: Fixed beta, varying dimensionless g
%plot the fixed points in the order of:
%beta = 72 degree, dimensionless g = [0.12, 0.2, 0.28, 0.36, 0.44]
h2 = figure();
hold on

colorData = 'maxAbsEigenValue';

plotStableFixedPointsPendulumSLIPPER('fixedPointData_Varing_none_082218_1529', colorData);
plotStableFixedPointsPendulumSLIPPER('fixedPointData_Varing_none_082218_1604', colorData);
plotStableFixedPointsPendulumSLIPPER('fixedPointData_Varing_none_082218_1617', colorData);
plotStableFixedPointsPendulumSLIPPER('fixedPointData_Varing_none_082218_1644', colorData);
plotStableFixedPointsPendulumSLIPPER('fixedPointData_Varing_none_082218_1649', colorData);

grid on

% font setup
set(gcf, 'DefaultTextFontName', 'Liberation Serif', 'DefaultTextFontSize', 18, ...
    'DefaultAxesFontName', 'Liberation Serif', 'DefaultAxesFontSize', 18, ...
    'DefaultLineLineWidth', 1, 'DefaultLineMarkerSize', 7.75)

% add parameter text
text(-0.0085,0.026,'$\beta = 72^o, \tilde{g} = [0.12, 0.2, 0.28, 0.36, 0.44]$','Interpreter','latex','FontSize',12)

% figure labels
xlabel('$\phi^*$', 'Interpreter', 'latex')
ylabel('$\dot{\phi^*}$', 'Interpreter', 'latex')

% colorbar setup
cbh = colorbar();
titleString = 'max(abs($\lambda$))';
ylabel(cbh, titleString, 'Interpreter', 'latex')

%% Delta-dimensiolness k plot: Fixed dimensionless g, varying beta
% plot the fixed points in the order of:
% dimensionless g = 0.2, beta = [66 68 70 72 74] (degree)
h3 = figure();
hold on

colorData = 'maxAbsEigenValue';
scaleByTotalMass = true;

plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_082218_2113', colorData, scaleByTotalMass);
plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_082218_2120', colorData, scaleByTotalMass);
plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_082218_2125', colorData, scaleByTotalMass);
plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_082218_2131', colorData, scaleByTotalMass);
plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_082218_2135', colorData, scaleByTotalMass);

plotStableFixedPointsSLIP('fixedPointData_Varing_beta_dutyFactorChecked_082118_1901');

axis([5, 20, 0, 0.3])
daspect([30, 1, 1])
grid on

% font setup
set(gcf, 'DefaultTextFontName', 'Liberation Serif', 'DefaultTextFontSize', 18, ...
    'DefaultAxesFontName', 'Liberation Serif', 'DefaultAxesFontSize', 18, ...
    'DefaultLineLineWidth', 1, 'DefaultLineMarkerSize', 7.75)

% add parameter text
text(8,0.32,'$\tilde{g} = 0.2, \beta = [66^o, 68^o, 70^o, 72^o, 74^o]$','Interpreter','latex','FontSize',12)

% figure labels
xlabel('$\tilde{k}$', 'Interpreter', 'latex')
ylabel('$\delta^*$', 'Interpreter', 'latex')

% colorbar setup
cbh = colorbar();
titleString = 'max(abs($\lambda$))';
ylabel(cbh, titleString, 'Interpreter', 'latex')

%% Phi-phid plot: Fixed beta, varying dimensionless g
%plot the fixed points in the order of:
%beta = 72 degree, dimensionless g = [0.12, 0.2, 0.28, 0.36, 0.44]
h2 = figure();
hold on

colorData = 'maxAbsEigenValue';

plotStableFixedPointsPendulumSLIPPER('fixedPointData_Varing_none_082218_2113', colorData);
plotStableFixedPointsPendulumSLIPPER('fixedPointData_Varing_none_082218_2120', colorData);
plotStableFixedPointsPendulumSLIPPER('fixedPointData_Varing_none_082218_2125', colorData);
plotStableFixedPointsPendulumSLIPPER('fixedPointData_Varing_none_082218_2131', colorData);
plotStableFixedPointsPendulumSLIPPER('fixedPointData_Varing_none_082218_2135', colorData);

grid on

% font setup
set(gcf, 'DefaultTextFontName', 'Liberation Serif', 'DefaultTextFontSize', 18, ...
    'DefaultAxesFontName', 'Liberation Serif', 'DefaultAxesFontSize', 18, ...
    'DefaultLineLineWidth', 1, 'DefaultLineMarkerSize', 7.75)

% add parameter text
text(-0.0085,0.0165,'$\tilde{g} = 0.2, \beta = [66^o, 68^o, 70^o, 72^o, 74^o]$','Interpreter','latex','FontSize',12)

% figure labels
xlabel('$\phi^*$', 'Interpreter', 'latex')
ylabel('$\dot{\phi^*}$', 'Interpreter', 'latex')

% colorbar setup
cbh = colorbar();
titleString = 'max(abs($\lambda$))';
ylabel(cbh, titleString, 'Interpreter', 'latex')
