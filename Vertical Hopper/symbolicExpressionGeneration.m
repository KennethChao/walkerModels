% script to generate closed-form Poincare map of Pitch dyanmics of a
% vertical hopper, with the combination of:
% 1) P control or PD control
% 2) Discrete or Continuous control during stance phase
% 
% Note for the Discrete control, the shifted horizontal position only
% updated at every touch down, then it will keep as constant for the whole
% stance phase.
%
% To see the expression or get its latex expression, use "Ctrl+Enter" to 
% run the section you want to choose.

%% Poincare Map of Discrete Pitch Angle P Control
clc;
clear;
syms a K T

M = [1-a^2*K T; 
    -2*a*K 1];


%% Poincare Map of Discrete Pitch Angle P Control

clc;
clear;
syms a K T C

M = [1-a^2*K T-a^2*C; 
    -2*a*K 1-a^2*C];

% Check for Jury Test of the Poincare stability of Discrete Pitch Angle PD Control

% detM = simplify(det(M));
% 
% pretty(detM)
% 
% trM = simplify(trace(M));
% 
% c2 = simplify(1-trM+detM);
% 
% pretty(c2)
% 
% c3 = simplify(1+trM+detM);
% 
% pretty(c3)

%% Poincare Map of Continuous Pitch Angle Control (P control)
clc;
clear;
syms a K T

A = [0 1; -2*K 0];
M = (expm(A*a)*[1 (T-a); 0 1])

pretty(M)
latex(M)

%% Poincare Map of Continuous Pitch Angle Control (PD control)
clc;
clear;
syms a K T C

A = [0 1; -2*K -2*C];
Mpd = (expm(A*a)*[1 (T-a); 0 1])

pretty(Mpd)
latex(Mpd)