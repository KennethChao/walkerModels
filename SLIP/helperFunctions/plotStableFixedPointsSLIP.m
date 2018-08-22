function plotStableFixedPointsSLIP(filename)
% PLOTSTABLEFIXEDPOINTSLIP helper function to plot stable fixed-point set of SLIP
% (for comparison between SLIP and SLIPPER)
%

% Load data
data = load(filename);
optParms = data.optParms;
disp(optParms);
result = data.result;

% Plot stable fixed points as 2D scatter plot (connecting with gray dash line)
    for k = 1:optParms.searchingVarLength 
        for i = 1: size(result.stableSolution,1)
            yData = result.stableSolution(i,:,k);
            colorData = result.stableData(i,:,k);
            plot(optParms.kVec, yData,'k--','linewidth',0.001,'Color',[0.8,0.8,0.8]);
            scatter(optParms.kVec, yData, 12,colorData,'filled');
            hold on
        end
    end

end