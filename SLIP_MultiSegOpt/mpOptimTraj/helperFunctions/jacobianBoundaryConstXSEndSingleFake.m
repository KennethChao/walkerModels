function jacobianStanceEnd = jacobianBoundaryConstXSEndSingleFake(l0,theta0,l0d,theta0d,lEnd,thetaEnd,lEndd,thetaEndd,g,k,beta)
%JACOBIANBOUNDARYCONSTXSENDSINGLEFAKE
%    JACOBIANSTANCEEND = JACOBIANBOUNDARYCONSTXSENDSINGLEFAKE(L0,THETA0,L0D,THETA0D,LEND,THETAEND,LENDD,THETAENDD,G,K,BETA)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    08-Sep-2018 11:26:37

t2 = sin(thetaEnd);
t3 = cos(thetaEnd);
t4 = cos(theta0);
t5 = l0d.*t4;
t6 = sin(theta0);
t7 = lEnd.*t2.*thetaEndd.*5.0;
t9 = lEndd.*t3.*5.0;
t10 = l0.*t6.*theta0d;
t8 = t5+t7-t9-t10;
t11 = lEndd.*t2.*5.0;
t12 = lEnd.*t3.*thetaEndd.*5.0;
t13 = l0d.*t6;
t14 = l0.*t4.*theta0d;
t15 = t11+t12+t13+t14;
jacobianStanceEnd = reshape([lEnd.*2.0-2.0,0.0,t2.*t8.*thetaEndd.*1.0e1,t3.*t15.*thetaEndd.*1.0e1,0.0,0.0,t8.*(t11+t12).*2.0,t15.*(t7-t9).*-2.0,0.0,0.0,t3.*t8.*-1.0e1,t2.*t15.*1.0e1,0.0,0.0,lEnd.*t2.*t8.*1.0e1,lEnd.*t3.*t15.*1.0e1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[4,6]);
