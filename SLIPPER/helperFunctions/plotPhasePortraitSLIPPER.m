function plotPhasePortraitSLIPPER(phasePortraitType, result, stableData, parms, plotParms)

%% Extract kinematic data from result
kinematicData = extractKinematicDatafromResult(result,parms);

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


grid on
hold on

%% Anime
% xPosVec = [xlegvec(1),xfvec(1),xvec(1)];
% yPosVec = [zlegvec(1),zfvec(1),zvec(1)];
%
% % xPosVec = [xfvec(1),xvec(1)];
% % yPosVec = [zfvec(1),zvec(1)];
%
% xPosLegVec = [xlegvec(1),xfvec(1)];
% yPosLegVec = [zlegvec(1),zfvec(1)];
%
% h0 = figure;
% h = plot(xPosVec,yPosVec,'-o');
% hold on
% h2 = plot(xPosVec,yPosVec,'r');
% axis equal
% axis([-0.5 3 0 1.2])
%
% % Parameters for animated gif
% filename = 'SLIPP_Animated.gif'; % Specify the output file name
% DT = 2*1e-2;
%
% for i = 2:4:length(zfvec)
%
%
%     xPosVec = [xlegvec(i),xfvec(i),xvec(i)];
%     yPosVec = [zlegvec(i),zfvec(i),zvec(i)];
% %     xPosVec = [xfvec(i),xvec(i)];
% %     yPosVec = [zfvec(i),zvec(i)];
%     xPosLegVec = [xlegvec(i),xfvec(i)];
%     yPosLegVec = [zlegvec(i),zfvec(i)];
%
%
%     set(h,'XData',xPosVec,'YData',yPosVec)
%     set(h2,'XData',xPosLegVec,'YData',yPosLegVec)
%
%     if (i>=length(xvecStance)+1)
%         set(h2,'Color','b')
%     end
%
%       frame = getframe(h0);
%       im = frame2im(frame);
%       [imind,cm] = rgb2ind(im,256);
%
%
%       if i == 2
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',DT);
%       else
%           imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',DT);
%       end
%       pause(0.1);
%
% end
end
