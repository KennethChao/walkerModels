function [position,isterminal,direction] = slipEventFcn(t,x)
position = x(1)-1; % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 1;   % The zero can be approached from either direction
end