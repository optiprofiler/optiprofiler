%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = dixon( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Dixon's tridiagonal quadratic  problem
%
%   Source: problem 156 (p. 51) in
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

pname   = 'dixon';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 2 )
         disp( [ ' ERROR in dixon: n = ', int2str(n),' should be > 1!' ] )
      end
   else
      n = 10;
   end
   varargout{1} = -ones( n, 1);                 % x0
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'QUR2-AY-V-0';                % class

case 'cpsstr'

   n          = varargin{1};
   eldom      = cell( n, 1 );
   eldom{ 1 } = [ 1 ];
   for iel = 2:n-1
      eldom{ iel } = [ iel-1 iel ];
   end
   eldom{ n } = [ n ];
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

   iel   = varargin{1};
   x     = varargin{2};
   n     = varargin{3};
   if ( iel == 1 || iel == n )
      riel = 1 - x(1);
      varargout{1} = riel^2;
      if ( nargout > 1 )
         varargout{2} =  - 2 * riel;
	 if ( nargout > 2 )
	    varargout{3} = 2;
	 end
      end
   else
      riel  = x(1) - x(2);
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = [ 1; -1 ];
         varargout{2} = 2 * Jiel * riel;
         if ( nargout > 2 )
            varargout{3} = 2 * Jiel * Jiel.';
         end
      end
   end

return

end