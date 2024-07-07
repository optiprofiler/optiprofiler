function [ f, g ] = CHARDIS0v( )
%
%  An attempt to rewrite the CHARDIS0 problem in Matlab
%
%  20 V 2024
%

%np1   = 3;
%x     =  [-5.000000e+00;  6.123234e-16;  2.500000e+00; -6.123234e-16; -0.000000e+00; 0.000000e+00]; % from pycutest
np1   = 20;
x     =  [ 4.72908621e+00;  1.62349735e+00;  3.73803399e+00;  2.90942864e+00; ...
           2.44687334e+00;  3.74521846e+00;  1.03362310e+00;  4.08168533e+00; ...
          -3.25971101e-01;  3.93388616e+00; -1.47993051e+00;  3.37390173e+00; ...
          -2.31701590e+00;  2.51695022e+00; -2.77728553e+00;  1.50299177e+00; ...
          -2.85525640e+00;  4.76458024e-01; -2.59568764e+00 ;-4.33143659e-01; ...
          -2.08296415e+00; -1.12724383e+00; -1.42585594e+00; -1.54889244e+00; ...
          -7.39965256e-01; -1.68695086e+00; -1.30388440e-01; -1.57355446e+00; ...
           3.23007220e-01; -1.27552667e+00;  5.75734903e-01; -8.81227872e-01; ...
           6.23005665e-01; -4.84904773e-01;  4.97798548e-01; -1.70894457e-01; ...
           2.63157895e-01; -6.44550947e-17;  0.00000000e+00;  0.00000000e+00 ];

n     = np1 + np1;
f     = 0;
g     = 0;
gscal = 0.01;  %  group scaling
for i = 1:np1
    for j = i+1:np1

       %  A new group O(I,J)
       
       fin = 0;
       gin = zeros( n, 1 );
       ii  = 2*i-1;
       jj  = 2*j-1;
       [ fij, gij ] = elt( x(ii), x(jj) );
       fin = fin + fij;
       gin(ii) = gin(ii) + gij(1);
       gin(jj) = gin(jj) + gij(2);
       ii = 2*i;
       jj = 2*j;
       [ fji, gji ] = elt( x(ii), x(jj) );
       fin = fin + fji;
       gin(ii) = gin(ii) + gji(1);
       gin(jj) = gin(jj) + gji(2);
       [ fa, ga ] = grp( fin );
       f = f + fa / gscal;
       g = g + ga * gin / gscal;
    end
end

return
end

function [ fel, gel ] = elt( x, y )
xmy    = x - y;
fel    = xmy*xmy;
gel(1) =  2*xmy;
gel(2) = -2*xmy;
return
end

function [ fa, ga ] = grp( alpha )
fa =  1/alpha;
ga = -1/(alpha*alpha);
return
end

