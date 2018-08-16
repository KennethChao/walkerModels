% clc; 
% close all

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
%%

% position and velocity of mf
%     zfvecStance = x(:,1) .* sin(x(:,3));
%     xfvecStance = -x(:,1) .* cos(x(:,3)); 
%     
%     xf0 = xfvecStance(end);
%     xfd0 = -xe(2) * cos(xe(3)) + xe(1) * xe(4) * sin(xe(3));
%     
%     % position and velocity of m
%     zvecStance = x(:,1) .* sin(x(:,3)) - rc*sin(x(:,5));
%     xvecStance = -x(:,1) .* cos(x(:,3)) + rc*cos(x(:,5)); 
%     
%     x0 = xvecStance(end);    
%     xd0 = -xe(2) * cos(xe(3)) + xe(1) * xe(4) * sin(xe(3)) - rc*sin(xe(5))*xe(6);
l = x(:,1)';
ld = x(:,2)';
theta = x(:,3)';
thetad = x(:,4)';
phi = x(:,5)';
phid = x(:,6)';

mbMotion = bodyMassCartesianMotionStance(l,theta,phi,ld,thetad,phid,rc);
mfMotion = frameMassCartesianMotionStance(l,theta,phi,ld,thetad,phid,rc);

xf = mfMotion(1,:);
zf = mfMotion(2,:);
xfd = mfMotion(3,:);
zfd = mfMotion(4,:);

xb = mbMotion(1,:);
zb = mbMotion(2,:);
xbd = mbMotion(3,:);
zbd = mbMotion(4,:);

comMotion = motionPointMass2COM(xf,zf,xb,zb,xfd,zfd,xbd,zbd,rc,mf);
%%
%     xc = xf0*mf/(mf+1) + x0/(mf+1);
%     xdc = xfd0*mf/(mf+1) + xd0/(mf+1);
% 
% 
% % position and velocity of mf
%     rc2f = rc*1/(1+mf);
%     rc2b = rc*mf/(1+mf);
% 
%     zfvecFlight = x2(2:end,1) + rc2f*sin(x2(2:end,3));
%     xfvecFlight = xc + xdc*t2(2:end) - rc2f*cos(x2(2:end,3)); 
    
    % position and velocity of m
%     zvecFlight = x2(2:end,1) - rc2b*sin(x2(2:end,3));
%     xvecFlight = xc + xdc*t2(2:end) + rc2b*cos(x2(2:end,3)); 

    
    zc = x2(:,1)';
    zcd = x2(:,2)';
    phi = x2(:,3)';
    phid = x2(:,4)';    
    xcd = ones(size(zc))*comMotion(3,end);
    xc = comMotion(1,end) + t2'.*xcd;
    
    
pointMassMotion = motionCOM2PointMass(xc,zc,phi,xcd,zcd,phid,rc,mf);    
    
%%
xfvecStance = xf;
xfvecFlight = pointMassMotion(1,2:end);

zfvecStance = zf;
zfvecFlight = pointMassMotion(2,2:end);


xfvec = [xfvecStance xfvecFlight];
zfvec = [zfvecStance zfvecFlight];
zfdvec = [zfd pointMassMotion(6,2:end)];
% xlegvec = [zeros(size(xfvecStance)); xfvecFlight + cos(beta)];
% zlegvec = [zeros(size(zfvecStance)); zfvecFlight - sin(beta)];
% 
xvecStance = xb;
xvecFlight = pointMassMotion(3,2:end);

zvecStance = zb;
zvecFlight = pointMassMotion(4,2:end);


xvec = [xvecStance xvecFlight];
zvec = [zvecStance zvecFlight];

xdvec = [xbd pointMassMotion(7,2:end)];

zdvec = [zbd pointMassMotion(8,2:end)];
% 
zcvec = [comMotion(2,:) zc(2:end)];

zcdvec = [comMotion(4,:) zcd(2:end)];
% 
tvec = [t;t2(2:end)+t(end)];
% 
% 
phivec = [x(:,5);x2(2:end,3)];
dphivec = [x(:,6);x2(2:end,4)];
% 
% vecCheck = [xfvecFlight-xvecFlight zfvecFlight-zvecFlight];
% for i = 1:size(vecCheck,1)
%     if norm(vecCheck(i,:))-rc>1e-8
%         error('pendulum length changed!')
%     end
% end


c = parula(50);

a = round((tvec(end) )/2.0*50 );
%% Figures

% plot([tvec; tvec],[xvec; yvec])
% 
% plot(tvec,xcvec)
% hold on
% plot(tvec,zcvec)
% % 
% plot(tvec,xfvec)
% hold on
% plot(tvec,zfvec)
% 
% plot(phivec,dphivec,'linewidth',2);%,

plot(zcvec,zcdvec,'linewidth',2,'color', c(a,:))
colorbar

% xlabel('position')
% ylabel('velocity')

xlabel('\phi')
ylabel('\phi dot')
zlabel('z')
% plot(zfvec,zfdvec,'linewidth',1.5,'color', c(a,:))
% plot(tvec,xvec)
% hold on
% plot(tvec,zvec)
% plot(tvec,xfvec)
% hold on
% plot(tvec,phivec)
% plot3(phivec,dphivec,zvec,'linewidth',1.5) %,'color', c(a,:)
% plot3(zvec,phivec,dphivec,'linewidth',2,'color', c(a,:))
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

