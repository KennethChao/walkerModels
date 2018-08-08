clc;
clear;
close all;
addpath('./SLIPPER')
addpath('./SLIP')
warning off;

%%

%ToDo
% Better result saving (Done)
% Saving file with a better naming and right location (Done)
% Save seperate plots (Done)
% Correct the printout msg ()
% Comparison to SLIP model

% Energy (optional)

%% Beta and dimensionless g
betaVec = 70 / 180 * pi;
% betaVec = (66:8:74) / 180 * pi;
% betaVec = 60:5:80;
% gVec = [0.025 0.05 0.1 0.21 0.46 0.66]
gVec = linspace(0.1, 0.4, 2);
% gVec = 0.25;

%% Fixed-point Finder Setup
optParms = fixedPointFinderOptions(gVec, betaVec);

%% Find fixed points of SLIPPER model
result = findFixedPointsSLIPPER(optParms);

%%
stringDateTime = datestr(now, 'mmddyy_HHMM');
cd ./SLIPPER/ret
fileName = sprintf('fixedPointData_Varing_%s_%s.mat', optParms.searchingVar, stringDateTime);

save(fileName, 'result', 'optParms');
cd ../../