%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = clplateb( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The clamped plate problem by Strang, Nocedal and Dax.
%   The problem comes from the discretization the following problem
%   in mechanics: a plate is clamped on one edge and loaded on the
%   opposite side.  The plate is the unit square.
%
%   In this version of the problem, the weight wght is distributed
%   equally along the upper edge, introducing a symmetry with respect
%   to the vertical axis.
%
%   The plate is clamped on its lower edge, by fixing the
%   corresponding variables to zero.
%
%   Source:
%   J. Nocedal,
%   "Solving large nonlinear systems of equations arising in mechanics",
%   Proceedings of the Cocoyoc Numerical Analysis Conference, Mexico,
%   pp. 132-141, 1981.
%
%   Ph. Toint 22 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'clplateb';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class ]

   if ( length( varargin ) )
      n     = varargin{1};
      nsq   = sqrt( n );
      error = ( n < 16 || nsq ~= round( nsq ) );
      if ( error )
         disp( [ ' ERROR in clplateb: n = ', int2str(n), ' but should be a square and at least 16!' ] )
      end
   else
      n   = 16;
      nsq = 4;
   end
%   nsqm1      = nsq - 1;
%   h          = 1 / nsqm1;
   x0         = zeros( n, 1 );
   xlower     = -Inf * ones( n, 1 );
   xupper     =  Inf * ones( n, 1 );
   xlower( 1:nsq ) = 0;
   xupper( 1:nsq ) = 0;
   varargout{1} = x0;
   varargout{2} = 'unknown';      % fstar
   varargout{3} = '';             % xtype
   varargout{4} = xlower;
   varargout{5} = xupper;
   varargout{6} = [];             % clower
   varargout{7} = [];             % cupper
   varargout{8} = 'OXR2-MN-V-0';  % class

case 'cpsstr'

   n     = varargin{1};
   nsq   = sqrt( n );
   nsqm1 = nsq - 1;
   eldom = cell( nsqm1^2+1, 1 );
   iel   = 1;
   for ii = 2:nsq
      for jj = 2:nsq
         kij = (ii-1)*nsq+jj;
         eldom{ iel } = [ kij kij-1 ];             % A(I,J)
	 iel = iel + 1;
	 eldom{ iel } = [ kij (ii-2)*nsq+jj ];     % B(I,J)
	 iel = iel + 1;
	 eldom{ iel } = [ kij kij-1 ];             % C(I,J)
	 iel = iel + 1;
	 eldom{ iel } = [ kij (ii-2)*nsq+jj ];             % D(I,J)
	 iel = iel + 1;
      end
   end
   eldom{ iel } = [ n-nsq+1:n ];                   % W
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
   nsq   = sqrt( n );
   nsqm1 = nsq - 1;
   wght  = -0.1;
   if ( iel == 4*nsqm1^2 + 1 )
      varargout{1} = (wght/nsqm1) * sum(x);
      if ( nargout > 1 )
         varargout{2} = (wght/nsqm1)*ones(size(x));
	 if ( nargout > 2 )
	    varargout{3} = 0;
	 end
      end
   else
      switch( mod(iel,4) )
      case { 1, 2 }                     % A(I,J) and B(I,J)
         riel = x(1)-x(2);
         varargout{1} = 0.5*riel^2;
	 if ( nargout > 1 )
	    Jiel = [ 1; -1 ];
	    varargout{2} =  Jiel * riel;
	    if ( nargout > 2 )
	       varargout{3} = Jiel*Jiel';
	    end
         end	    
      case { 3, 0 }                     % C(I,J) and D(I,J)
         riel = (x(1)-x(2))^2;
	 varargout{1} = 0.5*n*riel^2;
         if ( nargout > 1 )
	    Jiel = 2*(x(1)-x(2))*[ 1; -1 ];
	    varargout{2} = n*Jiel*riel;
	    if ( nargout > 2 )
	       Hiel = 2 * [ 1 -1; -1  1 ];
	       varargout{3} = n*(Jiel*Jiel'+riel*Hiel);
	    end
	 end      
      end
   end
   
end

return

end