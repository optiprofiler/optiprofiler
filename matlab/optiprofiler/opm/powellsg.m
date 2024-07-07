%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = powellsg( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Variable dimension  problem
%
%   Source: problem 13 in
%      J.J. More, B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%   Also problem 34 (p. 85) in
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   If the dimension is unspecified, the default n = 4 is chosen.
%
%   Ph. Toint 25 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'powellsg';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 4 || round( n/4 ) ~= n/4 )
         disp( [ ' ERROR in powellsg: n = ', int2str(n), ' must be > 3 and satisfy n = 4*m!' ] )
      end
   else
      n = 4;
   end
   for i = 1:4:n
      varargout{1}([i i+1 i+2 i+3 ]) = [ -3; -1; 0; 1 ];     % x0
   end
   varargout{1} = varargout{1}';
   varargout{2} = 0;                                         % fstar
   varargout{3} = '';                                        % xtype
   varargout{4} = [];                                        % xlower
   varargout{5} = [];                                        % xupper
   varargout{6} = [];                                        % clower
   varargout{7} = [];                                        % cupper
   varargout{8} = 'SUR2-AY-V-0';                             % class

case 'cpsstr'

   n     = varargin{1};
   eldom = cell( n, 1 );
   for iel = 1:n/4
      eldom{ 4*iel-3 } = [ 4*iel-3 4*iel-2 ];
      eldom{ 4*iel-2 } = [ 4*iel-1 4*iel   ];
      eldom{ 4*iel-1 } = [ 4*iel-2 4*iel-1 ];
      eldom{ 4*iel }   = [ 4*iel-3 4*iel   ];
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
   switch ( mod( iel, 4 ) )
   case 1
      riel = x(1) - 10 * x(2);
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = [ 1; -10 ];
	 varargout{2} = 2 * Jiel * riel;
	 if ( nargout > 2 )
	    varargout{3} = 2 * Jiel * Jiel.';
	 end
      end
   case 2
     riel = sqrt(5) * ( x(1) - x(2) );
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = sqrt(5) * [ 1; -1 ];
	 varargout{2} = 2 * Jiel * riel;
	 if ( nargout > 2 )
	    varargout{3} = 2 * Jiel * Jiel.';
	 end
      end
   case 3
      riel = ( x(1) - 2 * x(2) ) ^2;
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = 2 * ( x(1) - 2 * x(2) ) * [ 1; -2 ];
	 varargout{2} = 2 * Jiel * riel;
	 if ( nargout > 2 )
	    Hiel = 2 * [ 1, -2;
	                -2,  4 ];
	    varargout{3} = 2 * ( Jiel * Jiel.' + riel * Hiel );
	 end
      end
   case 0
     riel = sqrt(10) * ( x(1) - x(2) )^2;
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = 2 * sqrt(10) * ( x(1) - x(2) ) * [ 1; -1 ];
	 varargout{2} = 2 * Jiel * riel;
	 if ( nargout > 2 )
	    Hiel = 2 * sqrt(10) * [ 1, -1;
	                           -1,  1 ];
	    varargout{3} = 2 * ( Jiel * Jiel.' + riel * Hiel );
	 end
      end
   end

end

return

end