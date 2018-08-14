function dynamicsGEN(expr, write2File,normalizedMass,massScaling,normalizedLength)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

if nargin ==0
    write2File = true;
    normalizedMass = true;
    normalizedLength = true;
    massScaling = 'bodyMass';
    expr =  kinematicsGEN(true,massScaling);
end

syms l ld theta thetad phi phid
syms ldd thetadd phidd
syms rc mf mb k l0 g


zf = expr.zf;
xf = expr.xf;
zfd = expr.zfd;
xfd = expr.xfd;
zb = expr.zb;
xb = expr.xb;
zbd = expr.zbd;
xbd = expr.xbd;
%%
xVec = [l;theta;phi];
dxVec = [ld;thetad;phid];
ddxVec = [ldd;thetadd;phidd];

T = 1/2*mf*xfd^2 + 1/2*mf*zfd^2 + 1/2*mb*xbd^2 + 1/2*mb*zbd^2;

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
    variableVector = [l, theta, phi, ld, thetad, phid,g,k, mf, rc];
    matlabFunction(Mmat,'File','inertiaMatrix','Vars',variableVector);
    matlabFunction(bvec,'File','nonInertiaTerms','Vars',variableVector);
    cd ../
    end
    
end

end

