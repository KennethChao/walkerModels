function dx = dymModelDoubleSupportDimensionless(t, xVec, parms)
% Get parameters
g = parms.g;
k = parms.k;
leadingFootPos = parms.leadingFootPos;
trailingFootPos = parms.trailingFootPos;

% Get state variables
x = xVec(1);
xd = xVec(2);
y = xVec(3);
yd = xVec(4);
% Evalute dimensionless EOMs of the SLIP model
posDiffSpring1 = [x,y]-trailingFootPos;
posDiffSpring2 = [x,y]-leadingFootPos;

deltaSpring1 = norm(posDiffSpring1)^2-1;
deltaSpring2 = norm(posDiffSpring2)^2-1;

theta1 = atan2(posDiffSpring1(2),posDiffSpring1(1));
theta2 = atan2(posDiffSpring2(2),posDiffSpring2(1));


xdd = -deltaSpring1*k*cos(theta1)-deltaSpring2*k*cos(theta2);
ydd = -deltaSpring1*k*sin(theta1)-deltaSpring2*k*sin(theta2)-g;

dx = [xd; xdd; yd; ydd];

end