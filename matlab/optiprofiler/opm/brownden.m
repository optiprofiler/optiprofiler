%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = brownden( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The Brown-Dennis problem in 4 variables.  This is a nonlinear
%   least-squares problems with 20 groups.
%
%   Source: Problem 16 in
%      J.J. More', B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%   Alsso problem 30 in 
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 17 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'brownden';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 25; 5; -5; -1 ];            % x0
   varargout{2} = 85822.2;                      % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AN-4-0';                % class

case 'cpsstr'

   eldom = cell( 20, 1 );
   for iel = 1:20
      eldom{ iel } = [ 1 2 3 4 ];
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
   t    = iel / 5;
   sint = sin(t);
   ri1  = x(1) +  t * x(2) - exp(t);
   ri2  = x(3) + x(4)*sint - cos(t);
   riel = ri1^2 + ri2^2;
   varargout{1} = riel^2;
   if ( nargout > 1 )
      Jiel = 2 * [ ri1; t*ri1; ri2; sint*ri2 ];
      varargout{2} = 2 * Jiel * riel;
      if ( nargout > 2 )
         Hiel = 2 * [ 1, t  ,  0  ,  0    ;
	              t, t^2,  0  ,  0    ;
		      0, 0  ,  1  , sint  ;
		      0, 0  , sint, sint^2 ];
         varargout{3} = 2 * ( Jiel*Jiel.' + riel*Hiel );
      end
   end

end

return

end