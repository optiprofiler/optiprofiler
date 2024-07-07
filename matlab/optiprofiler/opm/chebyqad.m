%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = chebyqad( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The Chebychev quadrature problem in variable dimension, using the
%   exact formula for the shifted Chebyshev polynomials.  This is a
%   nonlinear least-squares problem with n groups. The Hessian is full.
%
%   Source: Problem 35 in
%      J.J. More', B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%   Also problem 58 in 
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   If the dimension is unspecified, the default n = 10 is chosen.
%
%   Ph. Toint 22 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'chebyqad';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 2 )
         disp( [ ' ERROR in chebyqad: n = ', int2str(n),' but should satisfy n > 1!' ] )
      end
   else
      n = 10;
   end
   varargout{1} = [ 1:n ]'/(n+1);               % x0
   varargout{2} = [ 0 0.002516873 0.00650395 0.0045729551 ]; % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AN-V-0';                % class

case 'cpsstr'

   n       = varargin{1};
   eldom   = cell( n, 1 );
   for iel = 1:n
      eldom{ iel } = [ 1:n ];
   end
   cpsstr.name  = pname;
   cpsstr.eldom = eldom;
   cpsstr.param = {};
   varargout{1} = cpsstr;
   
case 'objf'   % varargout = [ f, g, H ]

   x = varargin{1};
   if ( nargin > 2 && isstruct( varargin{2} ) && isfield( varargin{2}, 'eldom' ) );
      cpsstr = varargin{2};
   else
      cpsstr = problem( 'cpsstr', length( x ) );
   end
   varargout = opm_eval_cpsf( pname, 'elobjf', x, cpsstr.eldom, nargout, cpsstr.param{:} );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   iel  = varargin{1};
   x    = varargin{2};
   n    = length( x );
   if ( mod( iel, 2 ) )
      riel = 0;
   else
      riel = -1 / ( iel^2 - 1 );
   end
   if ( nargout > 1 )
      Jiel = zeros( n, 1 );
      if ( nargout > 2 )
         Hiel = zeros( n, n );
      end
   end
   for j = 1:n
      [ cj, gcj, Hcj ] = chebyqad( 'chebpol', iel, x(j) );
      riel = riel - cj / n;
      if ( nargout > 1 )
         Jiel( j ) = - gcj / n;
         if ( nargout > 2 )
	    Hiel( j, j ) = - Hcj / n;
	 end
      end
   end
   varargout{1} = riel^2;
   if ( nargout > 1 )
      varargout{2} = 2 * Jiel * riel;
      if ( nargout > 2 )
         varargout{3} = 2 * ( Jiel*Jiel.' + riel*Hiel );
      end
   end

case 'chebpol' % The i-th shifted Chebychev polynomial at x

   i  = varargin{1};
   x  = varargin{2};
   T0 = [ 1 2*x-1 ];
   U0 = [ 1 4*x-2 ];
   T1 = [ 0 2 ];
   T2 = [ 0 0 ];
   for j = 3:i+1
      n = j-1;
      t = 2*x-1;
      T0(j) = 2*t*T0(n) - T0(n-1);
      U0(j) = 2*t*U0(n) - U0(n-1);
      if ( nargout > 1 )
         T1(j) = 2*n*U0(j-1);
	 if ( nargout > 2 )
	    if ( t == 1 )
	       T2(j) = 4*( n^4 - n^2)/3;
	    elseif ( t == -1 )
	       T2(j) = 4*(-1)^n*((n^4 - n^2)/3);
	    else
	       T2(j) = 4*n*(j*T0(j)-U0(j))/(t^2-1);
	    end
         end
      end
   end
   varargout{1} = T0(end);%D
   if ( nargout > 1 )
      varargout{2} = T1(end);
      if ( nargout > 2 )
         varargout{3} = T2(end);
      end
   end
end

return

end