[ f, g ] = CHARDIS0v( )

np1 = 3
n   = np1 + np1;
x =  [-5.000000e+00;  6.123234e-16;  2.500000e+00; -6.123234e-16; -0.000000e+00; 0.000000e+00]

f = 0;
g = 0;
fin = 0;
gin = zeros( n, 1 );
for i = 1:np1
    ip1 = i+1
    for j = 1:np1
       [ fij, gij ] = elt( x(i), x(3+j) )
       fin = fin + fij;
       gin = gin + gij;
       [ fij, gij ] = elt( x(3+i), x(j) )
       fin = fin + fij;
       gin = gin + gij;
    end
    [ fa, ga ] = grp( fin )
    f = f + 100 * fa;
    g = g + 100* ga * gin;
end

return

end

function [ fel, gel ] = elt( x, y )
xmy    = x - y;
fel    = xmy*xmy:
gel(1) = 2*xmy:
gel(2) = -2*xmy;
return

end

function [ fa, ga ] = grp( fin )
fa = 1/fin;
ga = -1/(fin*fin);
return

end

