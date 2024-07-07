%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = bard( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Bard problem in 3 variables.
%   This function is a nonlinear least squares with 15 groups.  
%
%   Source: Problem 3 in
%      J.J. More', B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%   Also problem 16 in 
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 17 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'bard';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   n            = 3;
   varargout{1} = ones( n, 1 );                 % x0
   varargout{2} = [ 0.008215 17.4286 ];         % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AY-3-0';                % class

case 'cpsstr'

   eldom = cell( 15, 1 );
   for iel = 1:15
      eldom{ iel } = [ 1 2 3 ];
   end
   y = [ 0.14 0.18 0.22 0.25 0.29 0.32 0.35 0.39 0.37 0.58 0.73 0.16 1.34 2.10 4.39 ];
   cpsstr.name  = pname;
   cpsstr.eldom = eldom;
   cpsstr.param = { y };
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
   y    = varargin{3};
   u   = iel;
   v   = 16 - iel;
   w   = min( u, v );
   den = x(2)*v + x(3)*w;
   riel = x(1) + u / den - y( iel );
   varargout{1} = riel^2;
   if ( nargout > 1 )
      Jiel = [ 1; -u * v / den^2; -u * w / den^2 ];
      varargout{2} = 2 * Jiel * riel;
      if ( nargout > 2 )
	 Hiel = [ 0,               0,                    0       ;
	          0, 2 * u * v^2   / den^3, 2 * u * v * w / den^3;
	          0, 2 * u * v * w / den^3, 2 * u * w^2   / den^3 ];
         varargout{3} = 2 * ( Jiel * Jiel.' + riel * Hiel );
      end
   end
   
end

return

end