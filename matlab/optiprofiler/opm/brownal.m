%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = brownal( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The Brown almost linear problem in variable dimension.  This is a nonlinear
%   least-squares problems with n groups.
%
%   Source: Problem 27 in
%      J.J. More', B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%   Alsso problem 79 in 
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

pname   = 'brownal';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 1 )
         disp( [ ' ERROR in brownal: n = ', int2str(n),' but should satisfy n >= 1!' ] )
      end
   else
      n = 10;
   end
   if ( n < 2 )
      disp( [ ' ERROR in brownal: n = ',int2str(n), ' but should be > 1!' ] )
   end
   varargout{1} = 0.5 * ones( n, 1 );           % x0
   varargout{2} = [ 0 1 ];                      % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AN-V-0';                % class

case 'cpsstr'

   n     = varargin{1};
   eldom = cell( n, 1 );
   for iel = 1:n-1
      eldom{ iel } = [ 1:n ];
   end
   eldom{ n } = [ 1:n ];
   cpsstr.name  = pname;
   cpsstr.eldom = eldom;
   cpsstr.param = { n };
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
   n    = varargin{3};
   if ( iel < n )
      riel         = x(iel) + sum(x) - (n+1);
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel      = ones( n, 1);
	 Jiel(iel) = 2;
	 varargout{2} = 2 * Jiel* riel;
	 if ( nargout > 2 )
	    varargout{3} = 2 * Jiel * Jiel.';
	 end
      end
   else
      riel         = 1 - prod(x);
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = zeros( n, 1 );
         for j = 1:n
            imj     = setdiff( [ 1:n ], j );
            Jiel(j) = -prod( x(imj) );
	 end
	 varargout{2} = 2 * Jiel * riel;
	 if ( nargout > 2 )
	    for j = 1:n
               Hiel( j, j ) = 0;
	       imj          = setdiff( [ 1:n ], j );
               for jj = 1:j-1
                  imjmjj        = setdiff( imj, jj );
	          Hiel( jj, j ) = -prod( x(imjmjj) );
	          Hiel( j, jj ) = Hiel( jj, j );
	       end
	    end
	    varargout{3} = 2 * ( Jiel*Jiel.' + riel*Hiel );
	 end
      end
   end
end

return

end