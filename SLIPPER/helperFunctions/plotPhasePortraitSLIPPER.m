function plotPhasePortraitSLIPPER(phasePortraitType, result, stableData, parms, plotParms)

%% Extract kinematic data from result
kinematicData = extractKinematicDatafromResult(result, parms);

zfVec = kinematicData.zfVec;
zfdVec = kinematicData.zfdVec;

phiVec = kinematicData.phiVec;
phidVec = kinematicData.phidVec;

%% Color map
c = parula(plotParms.stableFixedPointNumber);

a = round((stableData - plotParms.colorMapMin)/(plotParms.colorMapMax - plotParms.colorMapMin)*plotParms.stableFixedPointNumber);
if a <= 0
    a = 1;
end

%% Plot figure
switch phasePortraitType
    
    case 'phasePotrait_zf'
        plot(zfVec, zfdVec, 'linewidth', 1.2, 'color', c(a, :));
        xlabel('z')
        ylabel('\dot z')
    case 'phasePotrait_phi'
        plot(phiVec, phidVec, 'linewidth', 1.2, 'color', c(a, :))
        xlabel('\phi')
        ylabel('\phi dot')
    case 'phasePotrait_phiVec&zf'
        plot3(phiVec, phidVec, zfVec, 'linewidth', 1.2, 'color', c(a, :))
        xlabel('\phi')
        ylabel('\phi dot')
        zlabel('z')
    case 'phasePotrait_zVec&phi'
        plot3(zfVec, zfdVec, phiVec, 'linewidth', 1.2, 'color', c(a, :)) %,
        xlabel('z')
        ylabel('\dot z')
        zlabel('phi')
end

end
