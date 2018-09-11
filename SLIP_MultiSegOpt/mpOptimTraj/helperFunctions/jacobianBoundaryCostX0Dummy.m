function jacobianStance0 = jacobianBoundaryCostX0Dummy(l0,theta0,l0d,theta0d,lEnd,thetaEnd,lEndd,thetaEndd,g,k,beta)
%JACOBIANBOUNDARYCOSTX0DUMMY
%    JACOBIANSTANCE0 = JACOBIANBOUNDARYCOSTX0DUMMY(L0,THETA0,L0D,THETA0D,LEND,THETAEND,LENDD,THETAENDD,G,K,BETA)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    10-Sep-2018 16:12:30

t2 = sin(theta0);
t3 = cos(theta0);
t4 = l0d.*t3;
t5 = cos(thetaEnd);
t6 = sin(thetaEnd);
t7 = lEnd.*t6.*thetaEndd.*5.0;
t9 = lEndd.*t5.*5.0;
t10 = l0.*t2.*theta0d;
t8 = t4+t7-t9-t10;
t11 = l0d.*t2;
t12 = l0.*t3.*theta0d;
t13 = lEndd.*t6.*5.0;
t14 = lEnd.*t5.*thetaEndd.*5.0;
t15 = t11+t12+t13+t14;
jacobianStance0 = reshape([t2.*t8.*theta0d.*-2.0,t3.*t15.*theta0d.*2.0,t8.*(t11+t12).*-2.0,t15.*(t4-t10).*2.0,t3.*t8.*2.0,t2.*t15.*2.0,l0.*t2.*t8.*-2.0,l0.*t3.*t15.*2.0,0.0,0.0,0.0,0.0],[2,6]);
