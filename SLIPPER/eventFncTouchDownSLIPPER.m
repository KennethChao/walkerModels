function [position, isterminal, direction] = touchDownEventFcnSLIPPER(t, x, parms)
% y = sin(beta)
beta = parms.beta;
mf = parms.mf;
rc = parms.rc;
rc2w = 1 / (1 + mf);
yw = x(1) + rc * rc2w * cos(x(3));

position = yw - sin(beta); % The value that we want to be zero
isterminal = 1; % Halt integration
direction = -1; % The zero can be approached from either direction
end