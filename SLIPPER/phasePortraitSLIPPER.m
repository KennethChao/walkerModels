function phasePortraitSLIPPER(phasePortraitType, ret, stableData, parms, plotParms)

%%
x = ret.x;
xe = ret.xe;
t = ret.t;
x2 = ret.x2;
xe2 = ret.xe2;
t2 = ret.t2;
rc = parms.rc;
mf = parms.mf;
beta = parms.beta;

%% Extract kinematic data

l = x(:, 1)';
ld = x(:, 2)';
theta = x(:, 3)';
thetad = x(:, 4)';
phi = x(:, 5)';
phid = x(:, 6)';

mbMotion = bodyMassCartesianMotionStance(l, theta, phi, ld, thetad, phid, rc);
mfMotion = frameMassCartesianMotionStance(l, theta, phi, ld, thetad, phid, rc);

xf = mfMotion(1, :);
zf = mfMotion(2, :);
xfd = mfMotion(3, :);
zfd = mfMotion(4, :);

xb = mbMotion(1, :);
zb = mbMotion(2, :);
xbd = mbMotion(3, :);
zbd = mbMotion(4, :);

comMotion = motionPointMass2COM(xf, zf, xb, zb, xfd, zfd, xbd, zbd, rc, mf);

zc = x2(:, 1)';
zcd = x2(:, 2)';
phi = x2(:, 3)';
phid = x2(:, 4)';
xcd = ones(size(zc)) * comMotion(3, end);
xc = comMotion(1, end) + t2' .* xcd;
pointMassMotion = motionCOM2PointMass(xc, zc, phi, xcd, zcd, phid, rc, mf);

%% Combine Stance and Flight kinematic data
xfvecStance = xf;
xfvecFlight = pointMassMotion(1, 2:end);

zfvecStance = zf;
zfvecFlight = pointMassMotion(2, 2:end);


xfvec = [xfvecStance, xfvecFlight];
zfvec = [zfvecStance, zfvecFlight];
zfdvec = [zfd, pointMassMotion(6, 2:end)];
%
xvecStance = xb;
xvecFlight = pointMassMotion(3, 2:end);

zvecStance = zb;
zvecFlight = pointMassMotion(4, 2:end);


xvec = [xvecStance, xvecFlight];
xdvec = [xbd, pointMassMotion(7, 2:end)];

zvec = [zvecStance, zvecFlight];
zdvec = [zbd, pointMassMotion(8, 2:end)];
%
zcvec = [comMotion(2, :), zc(2:end)];
zcdvec = [comMotion(4, :), zcd(2:end)];
%
tvec = [t; t2(2:end) + t(end)];
%
phivec = [x(:, 5); x2(2:end, 3)];
dphivec = [x(:, 6); x2(2:end, 4)];
%

%% Color map
c = parula(plotParms.stableFixedPointNumber);

a = round((stableData - plotParms.colorMapMin)/(plotParms.colorMapMax - plotParms.colorMapMin)*plotParms.stableFixedPointNumber);
if a <= 0
    a = 1;
end

%% PLot figure
switch phasePortraitType
    
    case 'phasePotrait_zf'
        
        plot(zfvec, zfdvec, 'linewidth', 1.2, 'color', c(a, :));
        xlabel('z')
        ylabel('\dot z')
    case 'phasePotrait_phi'
        plot(phivec, dphivec, 'linewidth', 1.2, 'color', c(a, :))
        xlabel('\phi')
        ylabel('\phi dot')
    case 'phasePotrait_phiVec&zf'
        
        plot3(phivec, dphivec, zvec, 'linewidth', 1.2, 'color', c(a, :))
        xlabel('\phi')
        ylabel('\phi dot')
        zlabel('z')
    case 'phasePotrait_zVec&phi'
        plot3(zcvec, zcdvec, phivec, 'linewidth', 1.2, 'color', c(a, :)) %,
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
