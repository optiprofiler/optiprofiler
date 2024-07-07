%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = brkmcc( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   A nonquadratic problem in 2 variables.
%
%   Source:  problem 85 in 
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 17 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'brkmcc';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 1; 2 ];                     % x0
   varargout{2} = [ -Inf, 0.16904 ];            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'OUR2-AN-2-0';                % class

case 'cpsstr'

   cpsstr.name  = pname;
   cpsstr.eldom = { [ 1 2 ] };
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

   x    = varargin{2};
   px   = -0.25*x(1)^2 - x(2)^2 + 1;
   hx   = x(1) - 2*x(2) + 1;
   varargout{1} = (x(1)-2)^2 + (x(2)-1)^2 + 1/( 25*px ) + 5*hx^2;
   if ( nargout > 1 )
      Jpx = [ -0.5*x(1); -2*x(2) ];
      Jhx = [     1    ;    -2   ];
      varargout{2}   = [ 2*(x(1)-2)-Jpx(1)/(25*px^2)+10*hx*Jhx(1);
                         2*(x(2)-1)-Jpx(2)/(25*px^2)+10*hx*Jhx(2)];
      if ( nargout > 2 )
         varargout{3}  = [ 2+1/(50*px^2)-Jpx(1)*x(1)/(25*px^3)+10,     -Jpx(2)*x(1)/(25*px^3)-20      ;
	                       -Jpx(2)*x(1)/(25*px^3)-20    , 2+2/(25*px^2)-4*Jpx(2)*x(2)/(25*px^3)+40 ];
      end
   end
end

return

end