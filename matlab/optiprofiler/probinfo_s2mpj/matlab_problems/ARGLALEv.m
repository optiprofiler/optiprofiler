function [ f, grad ] = ALGLALEv()

N = 4;
M = 6;
fact  = - 2.0/M;
fact2 = 1 - fact;
x =  ones(N,1);
g = zeros(M,1);
J = zeros(M,N);
for i = 1:N
   for j = 1:i-1
      g(i)   = g(i)   + fact*x(j);
      J(i,j) = J(i,j) + fact;
   end
   g(i)   = g(i)   + fact2*x(i);
   J(i,i) = J(i,i) + fact2;
   for j = i+1:N
      g(i)   = g(i)   + fact*x(j);
      J(i,j) = J(i,j) + fact;
   end
end

for i = N+1:M
   for j = 1:N
      g(i)   = g(i)   + fact*x(j);
      J(i,j) = J(i,j) + fact;
   end
end

ng = norm(g);
f    = 0.5*ng*ng;
grad = J'*g;

return

end

