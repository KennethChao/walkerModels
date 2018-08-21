function cost = costFun(t,x,parms)
% cost = costFun(u)
%
% Cost is the integral of torque-squared.
%     y = x(:, 1) .* sin(x(:, 3));

g = parms.g;
k = parms.k;

x0 = x(:,1);
xF = x(:,end);

dt = t(end)-t(end-1);

[velVec,deltaNew, deltaOld,costFly] = dymFlightDimensionless(x0,xF,dt,parms);

    y = -x(1, :) .* sin(x(3, :)); 
    yd = x(2, :) .* sin(x(3, :)) + x(1, :) .* x(4, :) .* cos(x(3, :));
%     x = -x(1, :) .* cos(x(:, 3)); %#ok<NASGU>
    xd = -x(2, :) .* cos(x(3, :)) + x(1, :) .* x(4, :) .* sin(x(3,: ));
cost = (yd.^2+xd.^2 )+ costFly/size(x,2) ;


end