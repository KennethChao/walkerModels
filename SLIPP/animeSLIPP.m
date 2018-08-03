clc; 
close all

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
    zfvecStance = x(:,1) .* sin(x(:,3));
    xfvecStance = -x(:,1) .* cos(x(:,3)); 
    
    xf0 = xfvecStance(end);
    xfd0 = -xe(2) * cos(xe(3)) + xe(1) * xe(4) * sin(xe(3));
    
    % position and velocity of m
    zvecStance = x(:,1) .* sin(x(:,3)) - rc*sin(x(:,5));
    xvecStance = -x(:,1) .* cos(x(:,3)) + rc*cos(x(:,5)); 
    
    x0 = xvecStance(end);    
    xd0 = -xe(2) * cos(xe(3)) + xe(1) * xe(4) * sin(xe(3)) - rc*sin(xe(5))*xe(6);


%%
    xc = xf0*mf/(mf+1) + x0/(mf+1);
    xdc = xfd0*mf/(mf+1) + xd0/(mf+1);


% position and velocity of mf
    rc2f = rc*1/(1+mf);
    rc2b = rc*mf/(1+mf);

    zfvecFlight = x2(2:end,1) + rc2f*sin(x2(2:end,3));
    xfvecFlight = xc + xdc*t2(2:end) - rc2f*cos(x2(2:end,3)); 
    
    % position and velocity of m
    zvecFlight = x2(2:end,1) - rc2b*sin(x2(2:end,3));
    xvecFlight = xc + xdc*t2(2:end) + rc2b*cos(x2(2:end,3)); 

%%

xfvec = [xfvecStance; xfvecFlight];
zfvec = [zfvecStance; zfvecFlight];

xlegvec = [zeros(size(xfvecStance)); xfvecFlight + cos(beta)];
zlegvec = [zeros(size(zfvecStance)); zfvecFlight - sin(beta)];

xvec = [xvecStance; xvecFlight];
zvec = [zvecStance; zvecFlight];

xcvec = xfvec*mf/(mf+1) + xvec/(mf+1);
zcvec = zfvec*mf/(mf+1) + zvec/(mf+1);


tvec = [t;t2(2:end)+t(end)];


phivec = [x(:,5);x2(2:end,3)];

vecCheck = [xfvecFlight-xvecFlight zfvecFlight-zvecFlight];
for i = 1:size(vecCheck,1)
    if norm(vecCheck(i,:))-rc>1e-8
        error('pendulum length changed!')
    end
end


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
% 
% plot(tvec,xvec)
% hold on
% plot(tvec,zvec)
% plot(tvec,xfvec)
% hold on
plot(tvec,phivec)

%% Anime
xPosVec = [xlegvec(1),xfvec(1),xvec(1)];
yPosVec = [zlegvec(1),zfvec(1),zvec(1)];

% xPosVec = [xfvec(1),xvec(1)];
% yPosVec = [zfvec(1),zvec(1)];

xPosLegVec = [xlegvec(1),xfvec(1)];
yPosLegVec = [zlegvec(1),zfvec(1)];

h0 = figure;
h = plot(xPosVec,yPosVec,'-o');
hold on
h2 = plot(xPosVec,yPosVec,'r');
axis equal
axis([-0.5 3 0 1.2])

% Parameters for animated gif
filename = 'SLIPP_Animated.gif'; % Specify the output file name
DT = 2*1e-2;

for i = 2:4:length(zfvec)
    
    
    xPosVec = [xlegvec(i),xfvec(i),xvec(i)];
    yPosVec = [zlegvec(i),zfvec(i),zvec(i)];
%     xPosVec = [xfvec(i),xvec(i)];
%     yPosVec = [zfvec(i),zvec(i)];    
    xPosLegVec = [xlegvec(i),xfvec(i)];
    yPosLegVec = [zlegvec(i),zfvec(i)];
    
    
    set(h,'XData',xPosVec,'YData',yPosVec)
    set(h2,'XData',xPosLegVec,'YData',yPosLegVec)
    
    if (i>=length(xvecStance)+1)
        set(h2,'Color','b')
    end
    
      frame = getframe(h0); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 

      
      if i == 2
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',DT); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',DT); 
      end    
      pause(0.1);
    
end

