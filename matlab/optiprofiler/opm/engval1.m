%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = engval1( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The ENGVAL1 problem.
%
%   Source: problem 31 in
%      Ph.L. Toint,
%      "Test problems for partially separable optimization and results
%      for the routine PSPMIN",
%      Report 83/4, Department of Mathematics, FUNDP (Namur, B), 1983.
%   Also problem 172 in
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   If the dimension is unspecified, the default n = 10 is chosen.
%
%   Ph. Toint and S. Gratton, 22 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'engval1';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 2 )
         disp( [ ' ERROR in engval1: n = ', int2str(n),' should be > 1!' ] )
      end
   else
      n = 10;
   end
   varargout{1} = 2*ones( n, 1);                % x0
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'OUR2-AY-V-0';                % class

case 'cpsstr'

   n     = varargin{1};
   eldom = cell( n-1, 1 );
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
   t   = x(1)^2 + x(2)^2;
   varargout{1} = t^2 - 4 * x(1) + 3;
   if ( nargout > 1 )
      Jt = 2 * [ x(1); x(2) ];
      varargout{2} = 2 * t * Jt - [ 4; 0 ];
      if ( nargout > 2 )
         varargout{3} = 2*Jt*Jt.' + 4*t*[ 1, 0; 0, 1];
      end
   end
   
return

end