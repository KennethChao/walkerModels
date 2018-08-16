function dx = dymModelStanceSLIPPendulum(t, x, parms)
% Get parameters
g = parms.g;
k = parms.k;
mf = parms.mf;
rc = parms.rc;
I = parms.I;
% Get state variables
l = x(1);
ld = x(2);
theta = x(3);
thetad = x(4);
phi = x(5);
phid = x(6);

Mmat = inertiaMatrix(l,theta,phi,ld,thetad,phid,g,k,mf,rc,I);
bvec = nonInertiaTerms(l,theta,phi,ld,thetad,phid,g,k,mf,rc,I);

xdd = Mmat\bvec;

ldd = xdd(1);
thetadd = xdd(2);
phidd = xdd(3);

dx = [ld; ldd; thetad; thetadd; phid; phidd];
end