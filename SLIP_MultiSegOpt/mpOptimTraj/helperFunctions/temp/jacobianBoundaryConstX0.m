function jacobianStance0 = jacobianBoundaryConstX0(l0,theta0,l0d,theta0d,lEnd,thetaEnd,lEndd,thetaEndd,g,k,beta)
%JACOBIANBOUNDARYCONSTX0
%    JACOBIANSTANCE0 = JACOBIANBOUNDARYCONSTX0(L0,THETA0,L0D,THETA0D,LEND,THETAEND,LENDD,THETAENDD,G,K,BETA)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    11-Sep-2018 11:17:49

t2 = l0d.^2;
t3 = l0.^2;
t4 = theta0d.^2;
t5 = t3.*t4;
t6 = t2+t5;
t7 = 1.0./t6;
t8 = t2+t5-1.0;
jacobianStance0 = reshape([l0d.*t7.*theta0d,l0.*t4.*t8.*4.0,1.0,0.0,-l0.*t7.*theta0d,l0d.*t8.*4.0,l0.*l0d.*t7,t3.*t8.*theta0d.*4.0,0.0,0.0,0.0,0.0],[2,6]);