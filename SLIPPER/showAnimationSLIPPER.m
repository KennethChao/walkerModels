% This is the script to show the animation for a (stable) fixed-point.
% Currently this code only support single beta and signle g (user need to 
% specify the k at line 18 for different beta or g values).
% 

addpath('./result')
addpath('./result/SLIPPER_PControl')
close all;
clc;
clear;

data = load('fixedPointData_Varing_none_081918_1313.mat');
optParms = data.optParms;
disp(optParms)
result = data.result;

% user input 1 (index for selecting one of varying g or varying beta)
k = 1;
if k > optParms.searchingVarLength
    error('index is exceeding the maximum value!');
end

% store mesh data with assigned 'k' index
resultBuf = nan(size(result.meshgridK));
for i = 1:optParms.sampledNumberK
    for j = 1:optParms.sampledNumberDelta
        resultBuf(j, i) = result.stableData(j, i, k).maxAbsEigenValue;
    end
end
stableIndices = find(~isnan(resultBuf));

% reshape mesh data to single vector
resultBuf = squeeze(reshape(resultBuf, 1, 1, []));
result.meshgridK = squeeze(reshape(result.meshgridK, 1, 1, []));
result.meshgridDelta = squeeze(reshape(result.meshgridDelta, 1, 1, []));
result.stablePhiReshape = [squeeze(reshape(result.stablePhi(1, :, :), 1, 1, [])), ...
    squeeze(reshape(result.stablePhi(2, :, :), 1, 1, [])), ...
    squeeze(reshape(result.stablePhi(3, :, :), 1, 1, []))];

% user input 2 to select the index of stable periodic motion
msg = sprintf('Choose a index for the SLIPPER motion (from 1 to %d)', length(stableIndices));
prompt = msg;
i = input(prompt);
if i > length(stableIndices)
    error('index is exceeding the maximum value!');
end

%%  Assign free varaibels and parameters
stableIndex = stableIndices(i);
parms = optParms;
parms.mode = 'simulationCheck';
parms.k = result.meshgridK(stableIndex);
parms.delta0 = result.meshgridDelta(stableIndex);

parms.beta = optParms.beta(k);
parms.g = optParms.g(k);

parms.controlMode = optParms.controlMode;
parms.controlGain = optParms.controlGain;

parms.optWeighting = optParms.optWeighting;

x = result.stablePhiReshape(stableIndex, :);
ret = oneStepSimulationSLIPPER(x, parms);

%% Extract kinematic data from result
kinematicData = extractKinematicDatafromResult(ret, parms);
xfVec = kinematicData.xfVec;
zfVec = kinematicData.zfVec;

xbVec = kinematicData.xbVec;
zbVec = kinematicData.zbVec;

xfootVec = kinematicData.xfootVec;
zfootVec = kinematicData.zfootVec;

%% Animetion
xPosVec = [xfootVec(1), xfVec(1), xbVec(1)];
yPosVec = [zfootVec(1), zfVec(1), zbVec(1)];

xPosLegVec = [xfootVec(1), xfVec(1)];
yPosLegVec = [zfootVec(1), zfVec(1)];

h0 = figure;
h = plot(xPosVec, yPosVec, '-o');
hold on
h2 = plot(xPosVec, yPosVec, 'r');
axis equal
axis([-0.5, 3, 0, 1.2])

% Parameters for animated gif
filename = 'SLIPPER_Animation.gif'; % Specify the output file name
DT = 2 * 1e-2;

for i = 2:4:length(zfVec)
    
    xPosVec = [xfootVec(i), xfVec(i), xbVec(i)];
    yPosVec = [zfootVec(i), zfVec(i), zbVec(i)];
    xPosLegVec = [xfootVec(i), xfVec(i)];
    yPosLegVec = [zfootVec(i), zfVec(i)];
    
    set(h, 'XData', xPosVec, 'YData', yPosVec)
    set(h2, 'XData', xPosLegVec, 'YData', yPosLegVec)
    
    % set the leg line color to blue in the flight phase
    if (i >= length(xfootVec) + 1)
        set(h2, 'Color', 'b')
    end
    
    frame = getframe(h0);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);

    % write image to gif file with assigned delayed time
    if i == 2
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', DT);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', DT);
    end
    pause(0.1);
    
end