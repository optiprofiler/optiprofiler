%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = tcontact( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The nonlinear contact problem.
%   This problem is a variant of the minimum surface problem.
%   Given a set of boundary conditions on the four sides of the square
%   and one flat obstacle in its middle, one must find the surface which
%   meets these boundary conditions and is of minimum area.
%
%   The unit square is discretized into squares, in turn discretized into 2 triangles.
%   The height funciton is linear on each triangle and is discretized using P1 
%   finite elements. 
%   Given these heights, the area above a square is
%   approximated by the sum of the area above each triangle 
%     S1(i,j) = 0.5*sqrt( 1 + ( a1(i,j)**2 + b1(i,j)^2 ) ) * h^2 
%   where
%     a1(i,j) = ( x( i  , j+1 ) - x(i, j   ) ) / h
%     b1(i,j) = ( x( i+1, j+1 ) - x(i, j+1 ) ) / h
%   and
%     S2(i,j) = 0.5*sqrt( 1 + ( a2(i,j)**2 + b2(i,j)^2 ) ) * h^2
%   where
%     a2(i,j) = ( x( i+1, j   ) - x( i  , j ) ) / h
%     b2(i,j) = ( x( i+1, j+1 ) - x( i+1, j ) ) / h
%
%   The starting point is the surface
%     x(i,j) = 1-(2*tx-1)^2.
%   where both
%     tx = (i-1)*h     and      ty = (j-1)*h.
%   In this problem, the boundary conditions are
%   given by
%     x(i,1) = 1 - (2*tx-1)^2
%     x(i,p) = 1 - (2*tx-1)^2
%   where
%     tx = (i-1)/(p-1)
%   and
%     x(1,j) = x(p,j) = 0
%   The obstacle is defined by
%     x(i,j) = 1
%   where both tx and ty are such that
%     | tx -1/2 | <= 1/4   and  | ty -1/2 | <= 1/4.
%   Since a simple projection of the surface defining the starting point
%   onto the feasible domain activates all variables in the obstacle region,
%   these are fixed in this version of the problem.
%
%   The number of variables must be a square and be at least 49.
%
%   Source:
%
%   Problem COPS 17 (p. 39) in
%   E. Dolan and J.J. Moré,
%   "Benchmarking Optimization Software with COPS",
%   techreport,  Argonne National Laboratory, ANL-MCS-246, 2000.
%
%   OPM input: Ph. Toint 9 VIII 2021.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'tcontact';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n     = varargin{1};
      nsq   = sqrt( n );
      error = ( n < 49 || nsq ~= round( nsq ) );
      if ( error )
         disp( [ ' ERROR in tcontact: n = ', int2str(n), ' but should be a square and at least 49!' ] )
      end
   else
      n   = 100;
      nsq = 10;
   end
   nsqm1  = nsq - 1;
   h      = 1 / nsqm1;
   x0     = zeros( n, 1 );
   xlower = -Inf * ones( n, 1 );
   xupper =  Inf * ones( n, 1 );
   for iy = 1:nsq
      if ( iy == 1 )
         for ix = 1:nsq
            tx        = (ix-1)*h;
            x0( ix ) = 1 - (2*tx-1)^2;
         end
      elseif ( iy == nsq )
         for ix = 1:nsq
            ii       = ix + (nsq-1)*nsq;
            t        = (ix-1)*h;
            x0( ii ) = 1 - (2*t-1)^2;
         end
      else
	 ii          = (iy-1)*nsq+1;
%        t           = (iy-1)*h;
         x0( ii )    = 0;
	 ii          = iy*nsq;
         x0( ii )    = 0;
      end
   end

   %  Fix the values on the obstacle define the starting point at free nodes.

   onobs = [];
   for i = 2:nsqm1
      for j = 2:nsqm1
 	 tx = (i-1)*h;
	 ty = (j-1)*h;
         ii = (j-1)*nsq+i;
	 if ( abs(tx-0.5) <= 0.25 &&  abs(ty-0.5) <= 0.25 )
	    x0( ii )     = 1;
	    onobs(end+1) = ii;
	 else
	    x0( ii ) = 1-(2*tx-1)^2;
	 end
      end
   end

   ind   = reshape([1:n],nsq,nsq);
   onbnd = unique([ind(1,1:nsq),ind(nsq,1:nsq),ind(1:nsq,nsq)',ind(1:nsq,1)']);
   xlower(onbnd) = x0(onbnd);
   xupper(onbnd) = x0(onbnd);
   xlower(onobs) = x0(onobs);
   xupper(onobs) = x0(onobs);
   varargout{1}  = x0;
   varargout{2}  = 2.480157;       % fstar for n = 128^2
   varargout{3}  = '';             % xtype
   varargout{4}  = xlower;
   varargout{5}  = xupper;
   varargout{6}  = [];             % clower
   varargout{7}  = [];             % cupper
   varargout{8}  = 'OXR2-MY-V-0';  % class

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

   x   = varargin{1};
   n   = length( x );
   nsq = sqrt( n );
   if ( nargin > 2 && isstruct( varargin{2} ) && isfield( varargin{2}, 'eldom' ) )
      cpsstr = varargin{2};
   else
      cpsstr = problem( 'cpsstr', length( x ) );
   end
   varargout = opm_eval_cpsf( pname, 'elobjf', x, cpsstr.eldom, nargout, n );
   
%   Zero derivatives of fixed variables.

   if ( nargout > 1 )
      ind        = reshape([1:n],nsq,nsq);
      onbnd      = unique([ind(1,1:nsq),ind(nsq,1:nsq),ind(1:nsq,nsq)',ind(1:nsq,1)']);
      varargout{2}(onbnd) = 0;
      if ( nargout > 1 )
         varargout{3}(onbnd,:) = 0;
         varargout{3}(:,onbnd) = 0;
      end
   end

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

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

   x = varargin{1};
   y = varargin{2};
   nsq   = sqrt(length(x));
   if ( round (nsq) ~= nsq )
      varargout{1} = NaN;
      disp( ' ERROR in tcontact-innerprod: dimension is not square!' )
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
%         lcg = 1/2*(c1+c2); % 1/2 because area of triangle is h^2/2
         lcg = (1/4)*(c1+c2); % 1/2 because area of triangle is h^2/2

         %  L2 part

         c1  = h2/24 * [ x1, t1, z1 ] * [ 2,1,1; 1,2,1; 1,1,2 ]*[ x2; t2; z2 ];
         c2  = h2/24 * [ x1, y1, z1 ] * [ 2,1,1; 1,2,1; 1,1,2 ]*[ x2; y2; z2 ];
%         lcf = c1+c2;
         lcf = (c1+c2)/2;

         %  Accumulate

         varargout{1} = varargout{1} + lcg + lcf;
      end
   end
   
case 'M-matrix'

   n   = varargin{1};
   nsq = sqrt(n);
   if ( round (nsq) ~= nsq )
      varargout{1} = NaN;
      disp( ' ERROR in tcontact-Mmatrix: dimension is not square!' )
      return
   end
   h2 = 1/(nsq-1)^2;
   varargout{1} = sparse( n, n);
   for i=1:nsq-1
      for j=1:nsq-1
         ix  = nsq*(i-1) + j;          it = ix + 1;
         iy  = ix + nsq;               iz = iy + 1;
         ind = [ ix, iy, iz, it];
         varargout{1}(ind,ind) = varargout{1}(ind,ind) ...
	       +   1/4  * [ 2 -1  0 -1 ; -1  2 -1  0 ; 0 -1  2 -1 ; -1  0 -1  2 ] ... % grad-grad part
               +  h2/48 * [ 4  1  2  1 ;  1  2  1  0 ; 2  1  4  1 ;  1  0  1  2 ];    % L2 part
%	       +   1/2  * [ 2 -1  0 -1 ; -1  2 -1  0 ; 0 -1  2 -1 ; -1  0 -1  2 ] ... % grad-grad part
%              +  h2/24 * [ 4  1  2  1 ;  1  2  1  0 ; 2  1  4  1 ;  1  0  1  2 ];    % L2 part
      end
   end

   %  Zero elements corresponding to boundary points

   ind   = reshape([1:n],nsq,nsq);
   onbnd = unique([ind(1,1:nsq),ind(nsq,1:nsq),ind(1:nsq,nsq)',ind(1:nsq,1)']);
   varargout{1}(onbnd,:) = 0;
   varargout{1}(:,onbnd) = 0;

end

return

end

