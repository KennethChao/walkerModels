function [position,isterminal,direction] = liftOffDoubleSupportEventFcn(t,x)
% l = 1 (normalized legnth)
l = sqrt(x(1)^2+x(3)^2);
position = l-1; % The value that we want to be zero
isterminal = 1;  % Halt integration 
direction = 1;   % The zero can be approached from either direction
end