    
close all

%%

% position and velocity of mf
    yfvecStance = x(:,1) .* sin(x(:,3));
%     yfd0 = xe(2) * sin(xe(3)) + xe(1) * xe(4) * cos(xe(3));
    xfvecStance = -x(:,1) .* cos(x(:,3)); 
    xfd0 = -xe(2) * cos(xe(3)) + xe(1) * xe(4) * sin(xe(3));
    
    % position and velocity of m
    yvecStance = x(:,1) .* sin(x(:,3)) - rc*sin(x(:,5));
%     yd0 = xe(2) * sin(xe(3)) + xe(1) * xe(4) * cos(xe(3)) - rc*cos(xe(5))*xe(6);
    xvecStance = -x(:,1) .* cos(x(:,3)) + rc*cos(x(:,5)); 
    xd0 = -xe(2) * cos(xe(3)) + xe(1) * xe(4) * sin(xe(3)) - rc*sin(xe(5))*xe(6);


%%
    xc = xf0*mf/(mf+1) + x0/(mf+1);
    xdc = xfd0*mf/(mf+1) + xd0/(mf+1);


% position and velocity of mf
    rc2f = 1/(1+mf);
    rc2b = mf/(1+mf);
%     yw = x(1) + rc2f*sin(x(3));

    yfvecFlight = x2(2:end,1) + rc*rc2f*sin(x2(2:end,3));
%     yfd0 = xe(2) * sin(xe(3)) + xe(1) * xe(4) * cos(xe(3));
    xfvecFlight = xc + xdc*t2(2:end) - rc*rc2f*cos(x2(2:end,3)); 
%     xfd0 = -xe(2) * cos(xe(3)) + xe(1) * xe(4) * sin(xe(3));
    
    % position and velocity of m
    yvecFlight = x2(2:end,1) - rc*rc2b*sin(x2(2:end,3));
%     yd0 = xe(2) * sin(xe(3)) + xe(1) * xe(4) * cos(xe(3)) - rc*cos(xe(5))*xe(6);
    xvecFlight = xc + xdc*t2(2:end) + rc*rc2b*cos(x2(2:end,3)); 
%     xd0 = -xe(2) * cos(xe(3)) + xe(1) * xe(4) * sin(xe(3)) - rc*sin(xe(5))*xe(6);




%%

xfvec = [xfvecStance; xfvecFlight];
yfvec = [yfvecStance; yfvecFlight];

xlegvec = [zeros(size(xfvecStance)); xfvecFlight + cos(beta)];
ylegvec = [zeros(size(yfvecStance)); yfvecFlight - sin(beta)];


xvec = [xvecStance; xvecFlight];
yvec = [yvecStance; yvecFlight];




xPosvec = [xlegvec(1),xfvec(1),xvec(1)];
yPosvec = [ylegvec(1),yfvec(1),yvec(1)];

h = plot(xPosvec,yPosvec,'-o');
axis([-0.5 5 0 1.2])
for i = 2:4:length(yfvec)
    
    
    xPosvec = [xlegvec(i),xfvec(i),xvec(i)];
    yPosvec = [ylegvec(i),yfvec(i),yvec(i)];
    
    set(h,'XData',xPosvec,'YData',yPosvec)
    pause(0.1)
    
end