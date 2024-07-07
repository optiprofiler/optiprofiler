%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = dixmaanf( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The Dixon-Maany problem.
%   The Hessian matrix is 7-diagonal with some widely separated.
%
%   The dimension n must satisfy n = 3*m for some positive integer m > 1.
%
%   Source:  
%      L.C.W. Dixon and Z. Maany,
%      "A family of test problems with sparse Hessians for unconstrained
%      optimization",
%      TR 206, Numerical Optimization Centre, Hatfield Polytechnic, 1988.
%   Also problem 18 in 
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   If the dimension is unspecified, the default n = 12 is chosen.
%
%   Ph. Toint 22 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'dixmaanf';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 3 || round( n/3 ) ~= n/3 )
         disp( [ ' ERROR in dixmaanf: n = ', int2str(n),' does not satisfy n = 3*m!' ] )
      end
   else
      n = 12;
   end
   varargout{1} = [ 2*ones( n, 1) ];            % x0
   varargout{2} = 1;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'OUR2-AY-V-0';                % class

case 'cpsstr'
   n     = varargin{1};
   m     = n / 3;
   eldom = cell( n, 1 );
   for i = 1:n
      irange = [ i ];
      if ( i < n )
         irange = [ irange, i+1 ];
      end
      if ( i <= 2*m )
         irange = [ irange, i+m ];
      end
      if ( i <= m )
         irange = [ irange, i+2*m ];
      end
      eldom{ i } = irange;
   end
   alpha = 1;  beta  = 0.625; gamma = 0.625; delta = 0.625; bigk  = [ 1, 0, 0, 1 ]; % case F
   cpsstr.name  = pname;
   cpsstr.eldom = eldom;
   cpsstr.param = { n, alpha, beta, gamma, delta, bigk };
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

   iel   = varargin{1};
   x     = varargin{2};
   n     = varargin{3};
   alpha = varargin{4};
   beta  = varargin{5};
   gamma = varargin{6};
   delta = varargin{7};
   bigk  = varargin{8};
   a     = alpha*(iel/n)^bigk(1);
   b     = beta *(iel/n)^bigk(2);
   c     = gamma*(iel/n)^bigk(3);
   d     = delta*(iel/n)^bigk(4);
   m     = n / 3;
   lx    = length( x );
   varargout{1}    = 0.5 * a * x(1)^2;
   if ( iel == 1 )
      varargout{1} = varargout{1} + 1;
   end
   if ( iel < n )
      varargout{1} = varargout{1} + b * x(1)^2 * ( x(2) + x(2)^2 )^2;
   end
   if ( iel <= 2*m )
      varargout{1} = varargout{1} + c * x(1)^2 * x(3)^4;
   end
   if ( iel <= m )
      varargout{1} = varargout{1} + d * x(1) * x(4);
   end
   if ( nargout > 1 )
      varargout{2} = [ a * x(1); zeros( lx-1, 1 ) ];
      if ( iel < n )
	 varargout{2}(1) = varargout{2}(1) + 2 * b * x(1) * ( x(2) + x(2)^2 )^2;
	 varargout{2}(2) = 2 *  b * x(1)^2 * ( x(2) + x(2)^2 ) * ( 1 + 2 * x(2) );
      end
      if ( iel <= 2*m )
  	 varargout{2}(1) = varargout{2}(1) + 2 * c * x(1) * x(3)^4;
	 varargout{2}(3) = 4 * c * x(1)^2 * x(3)^3;
      end
      if ( iel <= m )
	 varargout{2}(1) = varargout{2}(1) + d * x(4);
	 varargout{2}(4) = d * x(1);
      end
      if ( nargout > 2 )
         varargout{3} = zeros( lx, lx );
	 varargout{3}( 1, 1 ) = a;
	 if ( iel < n )
	    varargout{3}( 1, 1 ) = varargout{3}( 1, 1 ) + 2 * b * ( x(2) + x(2)^2 )^2;
	    varargout{3}( 1, 2 ) = 4 * b * x(1) * ( x(2) + x(2)^2 ) * ( 1 + 2 * x(2) );
	    varargout{3}( 2, 1 ) = varargout{3}( 1, 2 );
	    varargout{3}( 2, 2 ) = 2 * b * x(1)^2 * (1 + 6 * x(2) + 6 * x(2)^2 );
	 end
	 if ( iel <= 2*m )
            varargout{3}( 1, 1 ) = varargout{3}( 1, 1 ) + 2 * c * x(3)^4;
	    varargout{3}( 1, 3 ) = 8 * c * x(1) * x(3)^3;
	    varargout{3}( 3, 1 ) = varargout{3}( 1, 3 );
            varargout{3}( 3, 3 ) = 12 * c * x(1)^2 * x(3)^2;
	 end
	 if ( iel <= m )
	    varargout{3}( 1, 4 ) = d;
	    varargout{3}( 4, 1 ) = varargout{3}( 1, 4 );
	 end
      end
   end

end

return

end