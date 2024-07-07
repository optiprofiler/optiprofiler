%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = helix( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The extended "Helix" problem, featuring a multidimensional helical
%   valley.
%
%   Source: problem 7 in
%      J.J. More', B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%   Also problem 12 (p. 58) in
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   If the dimension is unspecified, the default n = 10 is chosen.
%
%   Ph. Toint and S. Gratton,  22 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'helix';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 3 )
         disp( [ ' ERROR in helix: n = ', int2str(n),' should be > 2!' ] )
      end
   else
      n = 10;
   end
   varargout{1} = [ -1; zeros( n-1, 1) ];         % x0
   varargout{2} = 0;                              % fstar
   varargout{3} = '';                             % xtype
   varargout{4} = [];                             % xlower
   varargout{5} = [];                             % xupper
   varargout{6} = [];                             % clower
   varargout{7} = [];                             % cupper
   varargout{8} = 'SUR2-AY-V-0';                  % class

case 'cpsstr'

   n     = varargin{1};
   eldom = cell( n-2, 1 );
   for iel = 1:n-2
      eldom{ iel } = [ 1 iel+1 iel+2 ];
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

   x = varargin{2};
   a = 1 / ( 2 * pi );
   r = sqrt( x(1)^2 + x(2)^2 );
   if ( x(1) > 0 )
      theta = a * atan( x(2)/x(1) );
   elseif ( x(1) < 0 )
      theta = 0.5 + a * atan( x(2)/x(1) );
   else
      theta = Inf;
   end
   r1 = 10 * ( x(3) - 10 * theta );
   r2 = 10 * ( r - 1 );
   r3 = x(3);
   varargout{1} = r1^2 + r2^2 + r3^2;
   if ( nargout > 1 )
      J1 = 10 * [ 10*a*x(2)/r^2;  -10*a*x(1)/r^2; 1 ];
      J2 = 10 * [     x(1)/r      ;         x(2)/r   ; 0 ];
      J3 =      [        0        ;            0     ; 1 ];
      varargout{2} = 2 * ( J1*r1 + J2*r2 + J3*r3 );
      if ( nargout > 2 )
% Revise this once the gradient is ok	 
	 H1 = 100 * [ -2*a*x(1)*x(2)/r^4   , a*(1-2*x(2)^2/r^2)/r^2, 0;
	             a*(1-2*x(2)^2/r^2)/r^2,  2*a*x(1)*x(2)/r^4    , 0;
	                      0            ,          0            , 0 ];
         H2 = 10 * [  (1-x(1)^2/r^2)/r     ,  -x(1)*x(2)/r^3       , 0;
	               -x(1)*x(2)/r^3      ,  (1-x(2)^2/r^2)/r     , 0;
	                      0            ,          0            , 0 ];
         varargout{3} = 2 * ( J1*J1.' + r1*H1 + J2*J2.' + r2*H2 + J3*J3.' );
      end
   end

return

end