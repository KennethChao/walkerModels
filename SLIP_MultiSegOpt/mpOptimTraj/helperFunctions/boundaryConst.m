function c = boundaryConst(l0,theta0,l0d,theta0d,lEnd,thetaEnd,lEndd,thetaEndd,g,k,beta)
%BOUNDARYCONST
%    C = BOUNDARYCONST(L0,THETA0,L0D,THETA0D,LEND,THETAEND,LENDD,THETAENDD,G,K,BETA)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    10-Sep-2018 10:54:12

t2 = l0d.^2+l0.^2.*theta0d.^2-1.0;
c = [angle(-(l0d+l0.*theta0d.*1i).*(cos(theta0)+sin(theta0).*1i)),t2.^2];
