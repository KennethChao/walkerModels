function dynamicsGEN(expr, write2File,normalizedMass,massScaling,normalizedLength)
%DYNAMICSGEN Derive and export dynamic equations of SLIPPER via Lagrange Mechanics
%   Input arguments
%
%   expr: output of kinematicsGEN, contains the motion of body and frame mass
%   normalizedMass: if true, the mass in the dynamics will be nomarlized
%   massScaling: the option for scaling masses, default is 'bodyMass'
%   normalizedMass: if true, the length in the dynamics will be nomarlized
%   by the leg length lo

if nargin ==0
    addpath('../../util')
    write2File = true; % export matlab function
    normalizedLength = true;
    normalizedMass = true; % normalized the mass with the given massScaling    
    massScaling = 'bodyMass';
    expr =  kinematicsGEN(write2File,massScaling);
end

% state variables
syms l ld theta thetad phi phid
syms ldd thetadd phidd
% system parameters
syms rc mf mb k l0 g I

% motion of frame mass
zf = expr.zf;
xf = expr.xf;
zfd = expr.zfd;
xfd = expr.xfd;
% motion of body mass
zb = expr.zb;
xb = expr.xb;
zbd = expr.zbd;
xbd = expr.xbd;

%% Lagrange mechanics
xVec = [l;theta;phi];
dxVec = [ld;thetad;phid];
ddxVec = [ldd;thetadd;phidd];

T = 1/2*mf*xfd^2 + 1/2*mf*zfd^2 + 1/2*mb*xbd^2 + 1/2*mb*zbd^2 + 1/2*I*phid^2;

V = 1/2*k*(l-l0)^2 + mf*g*(zf) + mb*g*(zb);

[Mmat,bvec] = LagrangeMechanics(T,V,xVec,dxVec,ddxVec);

if normalizedMass 
Mmat = massNomralization(Mmat,massScaling);
bvec = massNomralization(bvec,massScaling);
end

if normalizedLength
    Mmat = subs(Mmat,l0,1);
    bvec = subs(bvec,l0,1);
    
    if write2File
    cd ./autoGen
    variableVector = [l, theta, phi, ld, thetad, phid,g,k, mf, rc, I];
    matlabFunction(Mmat,'File','inertiaMatrix','Vars',variableVector);
    matlabFunction(bvec,'File','nonInertiaTerms','Vars',variableVector);
    cd ../
    end
    
end

end

