function [position, isterminal, direction] = touchDownWalkingEventFcn(t, x, beta)
% y = sin(beta)
position = x(1)*sin(x(3)) - sin(beta); % The value that we want to be zero
isterminal = 1; % Halt integration
direction = -1; % The zero can be approached from either direction
end