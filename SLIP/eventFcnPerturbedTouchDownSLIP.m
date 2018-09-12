function [position, isterminal, direction] = eventFcnPerturbedTouchDownSLIP(t, x, perturbation, beta)
% Event condition: y = sin(beta)
position = x(1) - sin(beta)+ perturbation; % The value that we want to be zero
isterminal = 1; % Halt integration
direction = -1; % Negative direction only
end