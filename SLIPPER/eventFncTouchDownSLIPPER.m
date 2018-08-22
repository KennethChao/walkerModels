function [position, isterminal, direction] = eventFncTouchDownSLIPPER(t, x, parms)
% Event condition: yf = sin(beta)
beta = parms.beta;
mf = parms.mf;
rc = parms.rc;
rc2f = 1 / (1 + mf);
yf = x(1) + rc * rc2f * cos(x(3));

position = yf - sin(beta); % The value that we want to be zero
isterminal = 1; % Halt integration
direction = -1; % Negative direction only
end