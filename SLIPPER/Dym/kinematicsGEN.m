

function expr =  kinematicsGEN()
%KINEMATICS Summary of this function goes here
%   Detailed explanation goes here

%%
    syms l(t) dl(t) theta(t) dtheta(t) phi(t) dphi(t)
    syms ddl(t) ddtheta(t) ddphi(t)

% initializeSymbolicVariables();
    syms rc mf

%% Kinematics

% mass of frame
zf = l(t)*sin(theta(t));
xf = -l(t)*cos(theta(t));

dzf = timeDiff(zf);
dxf = timeDiff(xf);

% body mass
zb = zf -rc*cos(phi);
xb = xf -rc*sin(phi) ;

dzb = timeDiff(zb);
dxb = timeDiff(xb);

expr.zf = zf;
expr.xf = xf;
expr.dzf = dzf;
expr.dxf = dxf;
expr.zb = zb;
expr.xb = xb;
expr.dzb = dzb;
expr.dxb = dxb;

%
% x = sym('x', [1 6]);
syms x
xf = variableSubstitution(xf)
zf = variableSubstitution(zf)

posStance = [xf;zf];
syms l dl theta dtheta phi dphi
matlabFunction(posStance,'File','myfile');

end

function dexp =  timeDiff(exp)    
syms l(t) dl(t) theta(t) dtheta(t) phi(t) dphi(t)
syms ddl(t) ddtheta(t) ddphi(t)

dexp = diff(exp,t);
dexp = subs(dexp,diff(l(t),t),dl(t));
dexp = subs(dexp,diff(theta(t),t),dtheta(t));
dexp = subs(dexp,diff(phi(t),t),dphi(t));

dexp = subs(dexp,diff(dl(t),t),ddl(t));
dexp = subs(dexp,diff(dtheta(t),t),ddtheta(t));
dexp = subs(dexp,diff(dphi(t),t),ddphi(t));
end


function newExpr =  variableSubstitution(expr)    
syms l(t) dl(t) theta(t) dtheta(t) phi(t) dphi(t)
syms ddl(t) ddtheta(t) ddphi(t)

syms l dl theta dtheta phi dphi

x = sym('x', [1 6]);
expr = subs(expr,l(t),l);
expr = subs(expr,dl(t),dl);
expr = subs(expr,theta(t),theta);
expr = subs(expr,dtheta(t),theta);
expr = subs(expr,phi(t),phi);
expr = subs(expr,dphi(t),dphi);

newExpr = expr;

end

