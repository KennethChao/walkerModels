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
clc;
%% Delta-dimensiolness k plot: Fixed beta, varying dimensionless g
% plot the fixed points in the order of:
% beta = 72 degree, dimensionless g = [0.12, 0.2, 0.28, 0.36, 0.44]
h1 = figure();
hold on

colorData = 'maxAbsEigenValue';
scaleByTotalMass = false;

% plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_082218_1529', colorData, scaleByTotalMass);
% plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_081918_1321', colorData, scaleByTotalMass);

%  plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_081918_1326', colorData, scaleByTotalMass);
%  plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_082218_1604', colorData, scaleByTotalMass);

%  plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_081918_1313', colorData, scaleByTotalMass);
%  plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_082218_1617', colorData, scaleByTotalMass);

%  plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_081918_1306', colorData, scaleByTotalMass);
%  plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_082218_1644', colorData, scaleByTotalMass);

%  plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_081918_1332', colorData, scaleByTotalMass);
%  plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_082218_1649', colorData, scaleByTotalMass);

% plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_082018_2202', colorData, scaleByTotalMass);
% plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_082218_2113', colorData, scaleByTotalMass);

% plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_081918_1321', colorData, scaleByTotalMass);
% plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_082218_2120', colorData, scaleByTotalMass);

% plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_082018_2212', colorData, scaleByTotalMass);
% plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_082218_2125', colorData, scaleByTotalMass);

% plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_082018_2218', colorData, scaleByTotalMass);
% plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_082218_2131', colorData, scaleByTotalMass);


% plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_082018_2222', colorData, scaleByTotalMass);
% plotStableFixedPointsSLIPPER('fixedPointData_Varing_none_082218_2135', colorData, scaleByTotalMass); 
 
 