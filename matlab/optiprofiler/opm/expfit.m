%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = expfit( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   A simple exponential fit in 2 variables
%
%   Source:
%   A.R. Conn, N. Gould and Ph.L. Toint,
%   "LANCELOT, a Fortran package for large-scale nonlinear optimization",
%   Springer Verlag, FUNDP, 1992.
%
%   SIF input: Ph. Toint, Jan 1991.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'expfit';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n ~= 2 )
         disp( [ ' ERROR in expfit: n = ', int2str(n), ' is not equal to 2!' ] )
      end
   else
      n = 2;
   end
   pones = ones(n/2,1);
   varargout{1}(1:2:n,1) =  pones;
   varargout{1}(2:2:n,1) = -pones;                % x0
   varargout{2} = 0;                              % fstar
   varargout{3} = '';                             % xtype
   varargout{4} = [];                             % xlower
   varargout{5} = [];                             % xupper
   varargout{6} = [];                             % clower
   varargout{7} = [];                             % cupper
   varargout{8} = 'SUR2-AN-2-0';                  % class

case 'cpsstr'

   eldom = cell( 10, 1 );
   for iel = 1:10  % number of points
      eldom{iel} = [ 1 2 ];
   end
   h = 0.25;              % stepsize
   cpsstr.name  = pname;
   cpsstr.eldom = eldom;
   cpsstr.param = { h };
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
   h     = varargin{3};
   switch( nargout )
   case 1
      riel = expfit( 'expf', x, iel*h );
      varargout{1} = riel^2;
   case 2
      [ riel, Jiel ] = expfit( 'expf', x, iel*h );
      varargout{1} = riel^2;
      varargout{2} = 2*Jiel*riel;
   case 3
      [ riel, Jiel, Hiel ] = expfit( 'expf', x, iel*h );
      varargout{1} = riel^2;
      varargout{2} = 2*Jiel*riel;
      varargout{3} = 2*(Jiel*Jiel' + riel*Hiel);
   end
   
case 'expf'

   x      = varargin{1};
   ih     = varargin{2};
   expxih = exp( x(2)*ih );
   varargout{1} = x(1)*expxih;
   if ( nargout > 1 )
      varargout{2} = [ expxih; x(1)*ih*expxih ];
      if ( nargout > 2 )
         varargout{3} = [      0         ih*expxih;
	                 ih*expxih   x(1)*ih^2*expxih ];
      end
   end

end

return

end
