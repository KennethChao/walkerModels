function newExpr =  massNomralization(expr,method)    

switch method
    case 'bodyMass'
        syms mb
        newExpr = subs(expr,mb,1);
    case 'totalMass'        
        syms mb mf
        newExpr = subs(expr,mb+mf,1);        
end