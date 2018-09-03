function jacobianStanceEnd = jacobianBoundaryConstXSEnd(lEnd,thetaEnd,lEndd,thetaEndd,g,k,beta)
%JACOBIANBOUNDARYCONSTXSEND
%    JACOBIANSTANCEEND = JACOBIANBOUNDARYCONSTXSEND(LEND,THETAEND,LENDD,THETAENDD,G,K,BETA)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    03-Sep-2018 14:32:28

t2 = sin(thetaEnd);
t3 = cos(thetaEnd);
jacobianStanceEnd = reshape([0.0,0.0,0.0,0.0,t3,-t2.*thetaEndd,-t2,-t3.*thetaEndd,0.0,0.0,0.0,0.0,0.0,-lEnd.*t2,-lEndd.*t2-lEnd.*t3.*thetaEndd,-lEnd.*t3,-lEndd.*t3+lEnd.*t2.*thetaEndd,0.0,0.0,0.0,0.0,0.0,0.0,t3,0.0,-t2,0.0,0.0,0.0,0.0,0.0,0.0,-lEnd.*t2,0.0,-lEnd.*t3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[9,6]);
