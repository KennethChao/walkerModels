function plotPhasePortraitsSLIPPER(phasePortraitType, result, stableData, parms, plotParms)
%PLOTPHASEPORTRAITSLIPPER helper function to plot phase portraits of stable fixed
% points
%   
%  Plotting options are availble for various 2D or 3D phase portraits.
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
        xlabel('$z_f$', 'Interpreter', 'latex')
        ylabel('$\dot{z}_f$', 'Interpreter', 'latex')
    case 'phasePotrait_phi'
        plot(phiVec, phidVec, 'linewidth', 1.2, 'color', c(a, :))
        xlabel('$\phi$', 'Interpreter', 'latex')
        ylabel('$\dot{\phi}$', 'Interpreter', 'latex')
    case 'phasePotrait_phiVec&zf'
        plot3(phiVec, phidVec, zfVec, 'linewidth', 1.2, 'color', c(a, :))
        xlabel('$\phi$', 'Interpreter', 'latex')
        ylabel('$\dot{\phi}$', 'Interpreter', 'latex')
        zlabel('$z_f$', 'Interpreter', 'latex')
    case 'phasePotrait_zVec&phi'
        plot3(zfVec, zfdVec, phiVec, 'linewidth', 1.2, 'color', c(a, :)) %,
        xlabel('$z_f$', 'Interpreter', 'latex')
        ylabel('$\dot{z}_f$', 'Interpreter', 'latex')
        zlabel('$\phi$', 'Interpreter', 'latex')
end

end
