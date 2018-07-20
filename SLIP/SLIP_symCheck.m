syms ld l t td phi phid rc

xb = ld*cos(t)-l*sin(t)*td-rc*sin(t+phi)*(td+phid);
yb = ld*sin(t)+l*cos(t)*td+rc*cos(t+phi)*(td+phid);

result = simplify(xb^2 + yb^2,'Criterion','preferReal');

pretty(result)

