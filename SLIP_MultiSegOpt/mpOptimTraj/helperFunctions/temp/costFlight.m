function costFlight = costFlight(xF,zF,xFd,zFd,h,g)
%COSTFLIGHT
%    COSTFLIGHT = COSTFLIGHT(XF,ZF,XFD,ZFD,H,G)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    11-Sep-2018 00:57:11

costFlight = h.*(g.*zF+xFd.^2+zFd.^2);
