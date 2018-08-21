function plotStableFixedPointsSLIPPER(filename, colorData, scaleByTotalMass)
%PLOTSTABLEFIXEDPOINTSSLIPPER helper function to plot stable fixed points 
%   Stable fixed points are ploted as a 2D contour plot.
data = load(filename);
optParms = data.optParms;
disp(optParms)
result = data.result;

for k = 1:optParms.searchingVarLength
    resultBuf = nan(size(result.meshgridK));
    for i = 1:optParms.sampledNumberK
        for j = 1:optParms.sampledNumberDelta
            switch colorData
                case 'maxAbsEigenValue'
                    resultBuf(j, i) = result.stableData(j, i, k).maxAbsEigenValue;
                case 'netWork'
                    resultBuf(j, i) = result.stableData(j, i, k).netWork;
                case 'dutyFactor'
                    resultBuf(j, i) = result.stableData(j, i, k).dutyFactor;
                case 'optimizedCost'
                    resultBuf(j, i) = result.stableData(j, i, k).fval;
                case 'dimensionlessStiffness'
                    resultBuf(j, i) = result.meshgridK(j, i);
                case 'dimensionlessGravityAccel'
                    resultBuf(j, i) = result.meshgridDelta(j, i);
                otherwise
                    error('unknown color data type');
            end
        end
    end
    
    if scaleByTotalMass
        bodyScale = 1 + optParms.mf;
    else
        bodyScale = 1;
    end
    
    contourf(result.meshgridK/bodyScale, result.meshgridDelta, resultBuf, 'LineStyle', ':');
    
end

end