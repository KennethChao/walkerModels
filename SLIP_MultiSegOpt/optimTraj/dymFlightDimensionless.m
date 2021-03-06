function [velVecF,deltaNew, deltaOld,costFly] = dymFlightDimensionless(xS,xF,dt,parms)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    % Convert states to Cartesian space
    
    g = parms.g;
    beta = parms.beta;
    
    yF0 = xF( 1) * sin(xF( 3));
    yFd0 = xF( 2) * sin(xF( 3)) + xF( 1) * xF( 4) * cos(xF( 3));
    xF0 = -xF( 1) * cos(xF( 3)); %#ok<NASGU>
    xFd0 = -xF( 2) * cos(xF( 3)) + xF( 1) * xF( 4) * sin(xF( 3));
    
    xF0 = [yF0, yFd0];
    
    yS0 = xS( 1) * sin(xS( 3));
    ySd0 = xS( 2) * sin(xS( 3)) + xS( 1) * xS( 4) * cos(xS( 3));
    xS0 = -xS( 1) * cos(xS( 3)); %#ok<NASGU>
    xSd0 = -xS( 2) * cos(xS( 3)) + xS( 1) * xS( 4) * sin(xS( 3));
    
    xS0 = [yS0, ySd0];
    
    
    % Solve ode
    if yFd0<0
        xFd0 = 1e1;
        yFdEnd = -1e1;

            te2 = 20;        
    else
%         [t, x2, te2, ae2, ie2] = ode45(dymFlight, tspan, xF0, optionsFlight); %#ok<ASGLU> % Runge-Kutta 4th/5th order ODE solver
        p = [-0.5*g +yFd0 yF0-sin(beta)];
        r = roots(p);
        
        if r(1)>0 && isreal(r(1))
            te2 = r(1);
        elseif r(2)>0 && isreal(r(2))
            te2 = r(2);
        else
            te2 = 50;
        end
        yFdEnd = yFd0 - g*te2;

    end
    velVecS = [xSd0, -ySd0];
    deltaOld = atan2(velVecS(2), velVecS(1));    
    
    % Get states of the next step
    velVecF = [xFd0, -yFdEnd];
    deltaNew = atan2(velVecF(2), velVecF(1));
    
    tVec = linspace(0,te2,round(te2/dt));
    xdFlight = ones(1,length(tVec))*xFd0;
    ydFlight = yFd0 - g*tVec ;
    
    
    costFly =sum(xdFlight.^2 + ydFlight.^2);
end

