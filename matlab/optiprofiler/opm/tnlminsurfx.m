%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = tnlminsurfx( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The nonlinear minimum surface problem.
%   This problem comes from the discretization of the minimum surface
%   problem on the unit square: given a set of boundary conditions on
%   the four sides of the square, one must find the surface which
%   meets these boundary conditions and is of minimum area.
%
%   The unit square is discretized into squares, in turn discretized into 2 triangles.
%   The height function is linear on each triangle and is discretized using P1 
%   finite elements. Given these heights, the area above a square is
%   approximated by the sum of the area above each triangle 
%     S1(i,j) = sqrt( 1 + ( a1(i,j)**2 + b1(i,j)^2 ) ) * (h^2/2) 
%   where
%     a1(i,j) = ( x( i  , j+1 ) - x(i, j   ) ) / h
%     b1(i,j) = ( x( i+1, j+1 ) - x(i, j+1 ) ) / h
%   and
%     S2(i,j) = sqrt( 1 + ( a2(i,j)**2 + b2(i,j)^2 ) ) * (h^2/2)
%   where
%     a2(i,j) = ( x( i+1, j   ) - x( i  , j ) ) / h
%     b2(i,j) = ( x( i+1, j+1 ) - x( i+1, j ) ) / h
%
%   In the Nonlinear Mininum Surface, the boundary conditions are
%   given by the following nonlinear functions:
%     x(i,1) = 1 + 8*t + 10*(1-t)^2
%     x(i,p) = 5 + 8*t + 10*(2-t)^2
%   where
%     t = (i-1)/(p-1)
%   and
%     x(1,j) = 1 + 4*t + 10*(1+t)^2
%     x(p,j) = 9 + 4*t + 10*t^2
%   where
%     t = (j-1)/(p-1).
%
%   The number of variables must be a square and be at least 9.
%
%   Source:
%
%   Ph. Toint, S. Gratton 28 VII 2021.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'tnlminsurfx';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n     = varargin{1};
      nsq   = sqrt( n );
      error = ( n < 9 || nsq ~= round( nsq ) );
      if ( error )
         disp( [ ' ERROR in tnlminsurfx: n = ', int2str(n), ' but should be a square and at least 9!' ] )
      end
   else
      n   = 16;
      nsq = 4;
   end
   nsqm1      = nsq - 1;
   h          = 1 / nsqm1;
   x0         = zeros( n, 1 );
   xlower     = -Inf * ones( n, 1 );
   xupper     =  Inf * ones( n, 1 );
   for iy = 1:nsq
      if ( iy == 1 )
         for ix = 1:nsq
            t        = (ix-1)*h;
            x0( ix ) = 1 + 8*t + 10*(1-t)^2;
         end
      elseif ( iy == nsq )
         for ix = 1:nsq
            ii       = ix + (nsq-1)*nsq;
            t        = (ix-1)*h;
            x0( ii ) = 5 + 8*t + 10*(2-t)^2;
         end
      else
	 ii          = (iy-1)*nsq+1;
	 t           = (iy-1)*h;
         x0( ii )    = 1 + 4*t + 10*(1+t)^2;
	 ii          = iy*nsq;
         x0( ii )    = 9 + 4*t + 10*t^2;
      end
   end
   ind   = reshape([1:n],nsq,nsq);
   onbnd = unique([ind(1,1:nsq),ind(nsq,1:nsq),ind(1:nsq,nsq)',ind(1:nsq,1)']);
   xlower(onbnd) = x0(onbnd);
   xupper(onbnd) = x0(onbnd);
   varargout{1} = x0;
   varargout{2} = '???';          % fstar
   varargout{3} = '';             % xtype
   varargout{4} = xlower;
   varargout{5} = xupper;
   varargout{6} = [];             % clower
   varargout{7} = [];             % cupper
   varargout{8} = 'OXR2-MY-V-0';  % class

case 'cpsstr'

   n     = varargin{1};
   nsq   = sqrt( n );
   nsqm1 = nsq - 1;
   eldom = cell( nsqm1^2, 1 );
   for ii = 1:nsqm1
      for jj = 1:nsqm1
         iel          = (ii-1)*nsqm1 + jj;
         ixl          = (ii-1)*nsq   + jj;
	 eldom{ iel } = [ ixl ixl+1 ixl+nsq ixl+nsq+1 ];
	 %                 x     t      y       z
      end
   end
   cpsstr.name  = pname;
   cpsstr.eldom = eldom;
   cpsstr.param = { n };
   varargout{1} = cpsstr;

case 'objf'   % varargout = [ f, g, H ]

   x = varargin{1};
   n = length( x );
   if ( nargin > 2 && isstruct( varargin{2} ) && isfield( varargin{2}, 'eldom' ) )
      cpsstr = varargin{2};
   else
      cpsstr = problem( 'cpsstr', n );
   end
   varargout = opm_eval_cpsf( pname, 'elobjf', x, cpsstr.eldom, nargout, n );
   
case 'elobjf' % varargout = [ fiel, giel, Hiel ]

%  iel   = varargin{1};
   x     = varargin{2};
   n     = varargin{3};
   nel   = (sqrt(n)-1)^2;
   r1    = 1 + nel*((x(2)-x(1))^2+(x(4)-x(2))^2);
   r2    = 1 + nel*((x(3)-x(1))^2+(x(4)-x(3))^2);
   S1    = sqrt(r1);
   S2    = sqrt(r2);
   varargout{1}  = 0.5*(S1 + S2)/nel;
   if ( nargout > 1 )
      p1  = 2*[ -(x(2)-x(1));x(2)-x(1)-(x(4)-x(2));0;x(4)-x(2)];
      p2  = 2*[ -(x(3)-x(1));0;x(3)-x(1)-(x(4)-x(3));x(4)-x(3)];
      varargout{2} = 0.5*((0.5/S1)*p1 + (0.5/S2)*p2);
      if ( nargout > 2 )
         H1 = 2*[ 1 -1 0 0; -1 2 0 -1;  0 0 0  0; 0 -1  0 1 ];
	 H2 = 2*[ 1 0 -1 0;  0 0 0  0; -1 0 2 -1; 0  0 -1 1 ];
         varargout{3} = 0.5*((0.5/S1)*H1-0.25*nel*p1*p1'/(r1*S1) + (0.5/S2)*H2-0.25*nel*p2*p2'/(r2*S2));
      end
   end

case 'innerprod'

   x   = varargin{1};
   y   = varargin{2};
   nsq = sqrt( length( x ) );
   if ( round (nsq) ~= nsq )
      varargout{1} = NaN;
      disp( ' ERROR in tnlminsurfx-innerprod: dimension is not square!' )
      return
   end
   nsqm1 = nsq - 1;
   h2    = 1/nsqm1^2;
   varargout{1} = 0;
   for i=1:nsqm1
      for j=1:nsqm1

	 y1 = x(i*nsq+j);      z1 = x(i*nsq+j+1);
         x1 = x((i-1)*nsq+j);  t1 = x((i-1)*nsq+j+1);
	 
	 y2 = y(i*nsq+j);      z2 = y(i*nsq+j+1);
         x2 = y((i-1)*nsq+j);  t2 = y((i-1)*nsq+j+1);
	 
         %  Grad-Grad part
	 
         c1  = [ t1-x1, z1-t1 ]*[ t2-x2; z2-t2 ];
         c2  = [ y1-x1, z1-y1 ]*[ y2-x2; z2-y2 ];
         lcg = (c1+c2)/2; 

         %  L2 part

         c1  = h2/24 * [ x1, t1, z1 ] * [ 2,1,1; 1,2,1; 1,1,2 ]*[ x2; t2; z2 ];
         c2  = h2/24 * [ x1, y1, z1 ] * [ 2,1,1; 1,2,1; 1,1,2 ]*[ x2; y2; z2 ];
         lcf = (c1+c2);

         %  Accumulate

         varargout{1} = varargout{1} + lcg + lcf;
      end
   end
   
case 'M-matrix'

   n   = varargin{1};
   nsq = sqrt(n);
   if ( round (nsq) ~= nsq )
      varargout{1} = NaN;
      disp( ' ERROR in tnlminsurfx-Mmatrix: dimension is not square!' )
      return
   end
   nsqm1 = nsq-1;
   h2    = 1/nsqm1^2;
   varargout{1} = sparse( n, n );
   for i=1:nsqm1
      for j=1:nsqm1
         ix  = nsq*(i-1) + j;          it = ix + 1;
         iy  = ix + nsq;               iz = iy + 1;
         ind = [ ix, iy, iz, it];
         varargout{1}(ind,ind) = varargout{1}(ind,ind) ...
	       +   1/2  * [ 2 -1  0 -1 ; -1  2 -1  0 ; 0 -1  2 -1 ; -1  0 -1  2 ] ... % grad-grad part
               +  h2/24 * [ 4  1  2  1 ;  1  2  1  0 ; 2  1  4  1 ;  1  0  1  2 ];    % L2 part
      end
   end

   %  Set the rows and columns corresponding to boundary points to the identity.

   ind   = reshape([1:n],nsq,nsq);
   onbnd = unique([ind(1,1:nsq),ind(nsq,1:nsq),ind(1:nsq,nsq)',ind(1:nsq,1)']);
   varargout{1}(onbnd,:) = 0;
   varargout{1}(:,onbnd) = 0;
   for k = 1:length(onbnd)
       varargout{1}(onbnd(k),onbnd(k)) = 1;
   end
end

return

end

