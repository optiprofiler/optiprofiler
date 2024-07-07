%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = morebv( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Moré's boundary value nonlinear least-squares problem.
%
%   Source: Problem 28 in
%      J.J. More', B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%   Also problem 83 (p. 75) in 
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   If the dimension is unspecified, the default n = 12 is chosen.
%
%   NOTE: the starting point specified in the above references has been
%         changed to [0;1, ..., 1; 0 ] to avoid early termination at x0 for 
%         large scale versions.
%
%   Ph. Toint 22 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'morebv';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 3 )
         disp( [ ' ERROR in morebv: n = ', int2str(n),' and should be > 2' ] )
      end
   else
      n = 12;
   end
%  t            = [1:n-2]'/(n-1);
%  varargout{1} = [ 0; t.*(t-1); 0 ];           % the standard x0
   varargout{1} = [ 0; ones(n-2,1);0 ];         % the modified x0
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [ 0; -Inf*ones( n-2, 1 ); 0 ];% xlower
   varargout{5} = [ 0;  Inf*ones( n-2, 1 ); 0 ];% xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SXR2-AY-V-0';                % class

case 'cpsstr'

   n = varargin{1};
   eldom   = cell( n-2, 1 );
   for iel = 1:n-2
      eldom{ iel } = [ iel iel+1 iel+2 ];
   end
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
   h    = 1 / ( n - 1 );
   t    = iel * h;
   riel = 2 * x(2) - x(1) - x(3) + 0.5*h^2*(x(2) + t + 1)^3;
   varargout{1} = riel^2;
   if ( nargout > 1 )
      Jiel = [ -1; 2 + (3/2)*h^2*(x(2) + t + 1)^2; -1 ];
      varargout{2} = 2 * Jiel * riel;
      if ( nargout > 2 )
	 Hiel = [ 0,           0         , 0;
	          0, 3*h^2*(x(2) + t + 1), 0;
	          0,           0         , 0];
	 varargout{3} = 2 * ( Jiel*Jiel.' + riel*Hiel );
      end
   end
   
end

return

end