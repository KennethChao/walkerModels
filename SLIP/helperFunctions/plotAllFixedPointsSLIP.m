function plotAllFixedPointsSLIP(filename,dataType,startingPointOption)
% PLOTALLFIXEDPOINTSSLIP helper function to plot all fixed-point set of SLIP
% (for comparison between SLIP and SLIPPER)
%

% Load data
data = load(filename);
optParms = data.optParms;
disp(optParms);
result = data.result;

for k = 1:optParms.searchingVarLength
    
    % Remove repeated fixed point solutions
    if strcmp(dataType,'eigenValue')
        trauncatedUnstableSolution = removeRepeatedFixedPointsSLIP(result.unstableData(:, :,k));
        trauncatedStableSolution = removeRepeatedFixedPointsSLIP(result.stableData(:, :,k));    
    elseif strcmp(dataType,'delta')
        trauncatedUnstableSolution = removeRepeatedFixedPointsSLIP(result.unstableSolution(:, :,k));
        trauncatedUnstableData = removeRepeatedFixedPointsSLIP(result.unstableData(:, :,k));    
        trauncatedStableSolution = removeRepeatedFixedPointsSLIP(result.stableSolution(:, :,k));        
        trauncatedStableData = removeRepeatedFixedPointsSLIP(result.stableData(:, :,k));    
    end

    hold on
    
    % Plot dots of fixed point solutions
    for i = 1: size(trauncatedUnstableSolution,1)
        plot(optParms.kVec, trauncatedUnstableSolution(i, :), 'ro','MarkerSize',7,'MarkerFaceColor' , [1 0 0] );
%         scatter(optParms.kVec,trauncatedUnstableSolution(i, :),30,trauncatedUnstableData(i, :),'filled')
    end
    for i = 1: size(trauncatedStableSolution,1)
%         plot(optParms.kVec, trauncatedStableSolution(i, :), 'bo','MarkerSize',5,'MarkerFaceColor' , [0 0 1] );
          scatter(optParms.kVec,trauncatedStableSolution(i, :),70,trauncatedStableData(i, :),'filled')
    end
    caxis([-1,1])
    % Sorting fixed points with the specified starting point option
    reshapedFixedPoints  = reshapeFixedPointsSLIP(trauncatedStableSolution,trauncatedUnstableSolution,optParms.kVec,optParms.deltaMax,startingPointOption);

    ax = gca;
    ax.ColorOrderIndex = 1;

    % Plotting line along the sorted fixed points
    plot(reshapedFixedPoints(2,:),reshapedFixedPoints(1,:))  
end

% Set figure axis limit
if strcmp(dataType,'eigenValue')
    axis([optParms.kMinFig, optParms.kMaxFig, -1.5, 2.5]) 
elseif strcmp(dataType,'delta')
    axis([optParms.kMinFig, optParms.kMaxFig, optParms.deltaMinFig, optParms.deltaMaxFig])
end