% This is the script to plot all fixed-point set of SLIP model
%
%   In the first three section, the results are used to compare to the
%   reference. The rest of sections are used to show how the stable fixed
%   points distributed in the fast running speed.
%
%   User-options:
%   ydataType: 'eigenValue' or 'delta'
%   startingPointOption: 'minDelta', 'maxNorm', or 'minNorm'
%   (option to determine the starting point of the line to connect fixed 
%   points)
%
%   Reference:
%   Shen Z, Seipel J. "A Piecewise-Linear Approximation of the Canonical 
%   Spring-Loaded Inverted Pendulum Model of Legged Locomotion.", ASME. J. 
%   Comput. Nonlinear Dynam., 2016.


clc;
close all;
clear;

addpath('./helperFunctions');
addpath('./result');
addpath('./result/Comparison2Paper');
addpath('./result/fastRunning');

%% Fixed dimensionless g = 0.46, varying beta = [66, 68, 70, 72, 74] degree
h1=   figure();

set(h1,'DefaultTextFontName','Liberation Serif','DefaultTextFontSize',18,...
'DefaultAxesFontName','Liberation Serif','DefaultAxesFontSize',18,...
'DefaultLineLineWidth',1,'DefaultLineMarkerSize',7.75)

ydataType = 'eigenValue';
startingPointOption = 'minDelta';

plotAllFixedPointsSLIP('fixedPointData_Varing_beta_082118_2247.mat', ydataType, startingPointOption);

xlabel('$\tilde{k}$','Interpreter','latex')
ylabel('$\lambda_2$','Interpreter','latex')
daspect([3, 1, 1])

%% Fixed dimensionless g = 0.46, varying beta = [66, 68, 70, 72, 74] degree
h2=   figure();

set(h2,'DefaultTextFontName','Liberation Serif','DefaultTextFontSize',18,...
'DefaultAxesFontName','Liberation Serif','DefaultAxesFontSize',18,...
'DefaultLineLineWidth',1,'DefaultLineMarkerSize',7.75)

ydataType = 'delta';
startingPointOption = 'minDelta';

plotAllFixedPointsSLIP('fixedPointData_Varing_beta_082118_2247.mat', ydataType, startingPointOption);

xlabel('$\tilde{k}$','Interpreter','latex')
ylabel('$\delta^*$','Interpreter','latex')
daspect([12, 1, 1])

%% Fixed beta = 72 degree, varying dimensionless g = [0.21, 0.46, 0.66, 0.86, 1.11, 1.31, 1.51]
h3=   figure();

set(h3,'DefaultTextFontName','Liberation Serif','DefaultTextFontSize',18,...
'DefaultAxesFontName','Liberation Serif','DefaultAxesFontSize',18,...
'DefaultLineLineWidth',1,'DefaultLineMarkerSize',7.75)

ydataType = 'delta';
startingPointOption = 'minDelta';

plotAllFixedPointsSLIP('fixedPointData_Varing_g_082118_2222.mat', ydataType, startingPointOption);

xlabel('$\tilde{k}$','Interpreter','latex')
ylabel('$\delta^*$','Interpreter','latex')
daspect([12, 1, 1])

%% Fixed dimensionless g = 0.1 (22.16 mph), varying beta = [66, 68, 70, 72, 74] degree
h4=   figure();

set(h4,'DefaultTextFontName','Liberation Serif','DefaultTextFontSize',18,...
'DefaultAxesFontName','Liberation Serif','DefaultAxesFontSize',18,...
'DefaultLineLineWidth',1,'DefaultLineMarkerSize',7.75)

ydataType = 'delta';
startingPointOption = 'maxNorm';

plotAllFixedPointsSLIP('fixedPointData_Varing_beta_082118_2349.mat', ydataType, startingPointOption);

xlabel('$\tilde{k}$','Interpreter','latex')
ylabel('$\delta^*$','Interpreter','latex')
daspect([15, 1, 1])

%% Fixed dimensionless g = 0.05 (31.33 mph), varying beta = [66, 68, 70, 72, 74] degree
h5=   figure();

set(h5,'DefaultTextFontName','Liberation Serif','DefaultTextFontSize',18,...
'DefaultAxesFontName','Liberation Serif','DefaultAxesFontSize',18,...
'DefaultLineLineWidth',1,'DefaultLineMarkerSize',7.75)

ydataType = 'delta';
startingPointOption = 'maxNorm';

plotAllFixedPointsSLIP('fixedPointData_Varing_beta_082118_2335.mat', ydataType, startingPointOption);

xlabel('$\tilde{k}$','Interpreter','latex')
ylabel('$\delta^*$','Interpreter','latex')
daspect([15, 1, 1])

%% Fixed dimensionless g = 0.025 (44.31 mph), varying beta = [66, 68, 70, 72, 74] degree
h6=   figure();

set(h6,'DefaultTextFontName','Liberation Serif','DefaultTextFontSize',18,...
'DefaultAxesFontName','Liberation Serif','DefaultAxesFontSize',18,...
'DefaultLineLineWidth',1,'DefaultLineMarkerSize',7.75)

ydataType = 'delta';
startingPointOption = 'maxNorm';

plotAllFixedPointsSLIP('fixedPointData_Varing_beta_082118_2340.mat', ydataType, startingPointOption);

xlabel('$\tilde{k}$','Interpreter','latex')
ylabel('$\delta^*$','Interpreter','latex')
daspect([15, 1, 1])