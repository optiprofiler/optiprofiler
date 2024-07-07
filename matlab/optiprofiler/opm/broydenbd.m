%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = broydenbd( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Broyden banded problem in variable dimension.  This is a nonlinear
%   least-squares problem with n groups.
%
%   Source: Problem 31 in
%      J.J. More', B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%   Alsso problem 73 in 
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

pname   = 'broydenbd';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 1 )
         disp( [ ' ERROR in broydenbd: n = ', int2str(n),' but should satisfy n >= 1!' ] )
      end
   else
      n = 10;
   end
   if ( n < 2 )
      disp( [ ' ERROR in broydenbd: n = ', int2str(n),' and should be > 1' ] )
   end
   varargout{1} = -ones( n, 1 );                % x0
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AY-V-0';                % class

case 'cpsstr'

   n       = varargin{1};
   eldom   = cell( n-2, 1 );
   for iel = 1:n
      eldom{ iel } = [ max(1,iel-5):min(n,iel+1) ];
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
   lx   = length( x );
   if ( iel <= 6 )
      ii = iel;
   else
      ii = 6;
   end
   riel = x(ii) * ( 2 + 5 * x(ii)^2 ) + 1;
   for j = 1:lx
      if ( ii ~= j )
         riel  = riel - x(j) * ( 1 + x(j) );
      end
   end
   varargout{1} = riel^2;
   if ( nargout > 1 )
      Jiel = zeros( lx, 1);
      for j = 1:lx
         if ( j == ii )
            Jiel(j) = 2 + 15 * x(j)^2;
         else
            Jiel(j) = - 1 - 2 * x(j);
	 end
      end
      varargout{2}  = 2 * Jiel * riel;
      if ( nargout > 2 )
         Hiel = zeros( lx, lx );
         for j = 1:lx
            if ( j == ii )
	       Hiel(j,j) = 30 * x(j);
	    else
	       Hiel(j,j) = -2;
	    end
	 end
         varargout{3} = 2*( Jiel*Jiel.' + riel*Hiel );
      end
   end
end

return

end