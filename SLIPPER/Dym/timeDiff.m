function dexpr =  timeDiff(expr,x,xd,xdd)    
gx = [x;xd];
dgx = [xd;xdd];

dim = length(x);
for i = 1:dim
    dexpr = jacobian(expr,gx)*dgx;
end

end