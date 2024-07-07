%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = argauss( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Gaussian problem in 3 variables.  This function is a nonlinear
%   least squares with 15 groups.
%
%   Source: problem 9 in
%      J.J. More, B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%   Also problem 28 in 
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 17 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'argauss';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 0.4; 1; 0 ];                % x0
   varargout{2} = 0.11279327696e-7;             % fstar
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
   y = [ 0.0009 0.0044 0.0175 0.0540 0.1295 0.2420 0.3521 0.3989 0.3521 0.2420 0.1295 0.0540 0.0175 0.0044 0.009 ];
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

   iel     = varargin{1};
   x       = varargin{2};
   y       = varargin{3};
   tielmx3 = 0.5 * ( 8 - iel ) - x(3);
   expiel  = exp( -0.5 * x(2) * tielmx3^2 );
   riel    = x(1) * expiel - y(iel);
   varargout{1} = riel^2;
   if ( nargout > 1 )
      Jiel = [ expiel; -0.5*x(1)*tielmx3^2*expiel; x(1)*x(2)*tielmx3*expiel ];
      varargout{2} = 2 *  Jiel * riel;
      if ( nargout > 2 )
         Hiel = expiel * ...
	    [          0    ,             -0.5*tielmx3^2           ,  x(2)*tielmx3;
	      -0.5*tielmx3^2,        0.25*x(1)*tielmx3^4           , (x(1)*tielmx3)*(1-0.5*x(2)*tielmx3^2);
	       x(2)*tielmx3 , (x(1)*tielmx3)*(1-0.5*x(2)*tielmx3^2), (x(1)*x(2))*(-1+x(2)*tielmx3^2)      ];
	    varargout{3} = 2 * ( Jiel * Jiel.' + riel * Hiel );
      end
   end

end

return

end