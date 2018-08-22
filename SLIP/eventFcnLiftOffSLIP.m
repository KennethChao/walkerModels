function [position,isterminal,direction] = eventFcnLiftOffSLIP(t,x)
% Event condition: l = 1 (normalized legnth)
position = x(1)-1; % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 1;   % Positive direction only
end