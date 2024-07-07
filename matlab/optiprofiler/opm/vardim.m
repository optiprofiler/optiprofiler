%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = vardim( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Variable dimension  problem
%
%   Source: problem 25 in
%      J.J. More, B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%   Also problem 72 (p. 98) in
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

pname   = 'vardim';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 2 )
         disp( [ ' ERROR in vardim: n = ', int2str(n), ' must be > 1!' ] )
      end
   else
      n = 10;
   end
   varargout{1} = ones( n, 1) - [ 1:n ]'/n;     % x0
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AY-V-0';                % class

case 'cpsstr'

   n     = varargin{1};
   eldom = cell( n+1, 1 );
   for iel = 1:n
      eldom{ iel } = [ iel ];
   end
   eldom{ n+1 } = [ 1:n ];
   eldom{ n+2 } = [ 1:n ];
   cpsstr.name  = pname;
   cpsstr.eldom = eldom;
   cpsstr.param = { n };
   varargout{1} = cpsstr;

case 'objf'   % varargout = [ f, g, H ]

   x = varargin{1};
   if ( nargin > 2 && isstruct( varargin{2} ) && isfield( varargin{2}, 'eldom' ) )
      cpsstr = varargin{2};
   else
      cpsstr = problem( 'cpsstr', length( x ) );
   end
   varargout = opm_eval_cpsf( pname, 'elobjf', x, cpsstr.eldom, nargout, cpsstr.param{:} );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   iel  = varargin{1};
   x    = varargin{2};
   n    = varargin{3};
   if ( iel <= n )
      riel = x(1)-1;
      varargout{1} = riel^2;
      if ( nargout > 1 )
	 varargout{2} = 2 * riel;
         if ( nargout > 2 )
	    varargout{3} = 2;
	 end
      end
   elseif ( iel < n+2 )
      riel = sum( [1:n]*(x-1) );
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = [ 1:n ]';
	 varargout{2} = 2 * Jiel * riel;
         if ( nargout > 2 )
	    varargout{ 3 } = 2 * ( Jiel*Jiel.' );
	 end
      end
   else
      t    = sum( [1:n]*(x-1) );
      riel = t^2;
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = 2 * t * [ 1:n ]';
	 varargout{2} = 2 * Jiel * riel;
         if ( nargout > 2 )
	    Hiel = 2 * [1:n].'*[1:n];
	    varargout{3} = 2 * ( Jiel*Jiel.' + riel*Hiel );
	 end
      end
   end
end

return

end