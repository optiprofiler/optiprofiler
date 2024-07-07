%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = edensch( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The extended Dennis-Schnabel problem, as defined by Li.
%
%   Source:
%      "The secant/finite difference algorithm for solving sparse
%      nonlinear systems of equations",
%      SIAM Journal on Optimization, (to appear), 1990.
%
%   If the dimension is unspecified, the default n = 10 is chosen.
%
%   Ph. Toint and S. Gratton, 22 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'edensch';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 1 )
         disp( [ ' ERROR in edensch: n = ', int2str(n),' should be > 0!' ] )
      end
   else
      n = 10;
   end
   varargout{1} = [ 8*ones( n, 1) ];            % x0
   varargout{2} = 'unknown';                    % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'OUR2-AY-V-0';                % class

case 'cpsstr'

   n   = varargin{1};
   for iel = 1:n-1
      eldom{ iel } = [ iel iel+1 ];
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

   x   = varargin{2};
   aux = x(1) * x(2) - 2 * x(2);
   varargout{1}  = ( x(1)- 2 )^4 + aux^2 + ( x(2) + 1 )^2;
   if ( nargout > 1 )
      varargout{2} = [4*(x(1)- 2)^3+2*aux*x(2); 2*aux*(x(1)-2)+2*(x(2)+1)];
      if ( nargout > 2 )
         varargout{3} = [12*(x(1)- 2)^2+2*x(2)^2, 2*aux+2*(x(1)-2)*x(2);
                          2*aux+2*(x(1)-2)*x(2)  , 2*(x(1)-2)*(x(1)-2)+2 ];
      end
   end

return

end