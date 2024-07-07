%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = bratu1d( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Bratu's problem in one dimension, according to Osborne.
%
%   Source: problem 121 (p. 99) in
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 27 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'bratu1d';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 3 )
         disp( [ ' ERROR in bratu1d: n = ', int2str(n), ' must be > 2!' ] )
      end
   else
      n = 30;
   end
   varargout{1} = (1/n)*ones( n , 1 );          % x0
   varargout{2} = 'unknown';                    % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [ 0; -Inf*ones(n-2,1); 0 ];   % xlower
   varargout{5} = [ 0;  Inf*ones(n-2,1); 0 ];   % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AY-V-0';                % class

case 'cpsstr'

   n = varargin{1};
   eldom = cell( 3*(n-2), 1 );
   for iel = 1:n-2
      eldom{ 3*iel-2 } = [ iel+1 ];
      eldom{ 3*iel-1 } = [ iel iel+1 ];
      eldom{ 3*iel   } = [ iel+1 iel+2 ];
   end
   cpsstr.name  = pname;
   cpsstr.eldom = eldom;
   varargout{1} = cpsstr;
   
case 'objf'   % varargout = [ f, g, H ]

   x = varargin{1};
   if ( nargin > 2 && isstruct( varargin{2} ) && isfield( varargin{2}, 'eldom' ) );
      cpsstr = varargin{2};
   else
      cpsstr = problem( 'cpsstr', length( x ) );
   end
   n = length( x );
   varargout = opm_eval_cpsf( pname, 'elobjf', x, cpsstr.eldom, nargout, n );

   %  Zero the derivatives of fixed variables.
   
   if ( nargout > 1 )
      varargout{2}( [ 1, n ] ) = 0;
      if ( nargout > 2 )
         varargout{3}( [ 1, n ], : ) = zeros( 2, n );
         varargout{3}( :, [ 1, n ] ) = zeros( n, 2 );
      end
   end

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   iel    = varargin{1};
   x      = varargin{2};
   n      = varargin{3};
   lambda = -3.4;
   h      = 1/(n-1);
   switch ( mod( iel, 3 ) );
   case 1
      varargout{1} = (2/h)*x(1)^2;
      if ( nargout > 1 )
         varargout{2} = (4/h)*x(1);
         if ( nargout > 2 )
	    varargout{3} = 4/h;
	 end
      end
   case 2
      varargout{1} = (2/h)*x(1)*x(2);
      if ( nargout > 1 )
         varargout{2} = (2/h)*[ x(2); x(1)];
         if ( nargout > 2 )
	    varargout{3} = (2/h)*[ 0, 1; 1, 0 ];
	 end
      end
   case 0
      e1  = exp(x(1));
      e2  = exp(x(2));
      dex = e2-e1;
      dx  = x(2)-x(1);
      varargout{1} = 2*lambda*h*dex/dx;
      if ( nargout > 1 )
         dx2    = dx^2;
	 dexdx2 = dex/dx2;
         varargout{2} = 2*lambda*h*[  dexdx2-e1/dx; -dexdx2+e2/dx ];
         if ( nargout > 2 )
	    dx3     = dx^3;
	    dexdx3  = dex/dx3;
	    H(1,1) = -e1/dx-2*(e1/dx2 -dexdx3);
	    H(2,1) = (e1+e2)/dx2-2*dexdx3;
	    H(1,2) = H(2,1);
	    H(2,2) = e2/dx - 2*(e2/dx2 - dexdx3 );
	    varargout{3} = 2*lambda*h*H;
	 end
      end
   end

end

return

end