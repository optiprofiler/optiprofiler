%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = biggs6( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Biggs problem in 6 variables.
%   This function is a nonlinear least squares with 13 groups.  
%
%   Source: Problem 18 in
%      J.J. More', B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%   Also problem 21 in 
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 17 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'biggs6';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 1; 2; 1; 1; 1; 1 ];         % x0
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AN-6-0';                % class

case 'cpsstr'

   eldom = cell( 13, 1 );
   for iel = 1:13
      eldom{ iel } = [ 1:6 ];
      ti = iel / 10;
      y( iel )     = exp( -ti ) - 5*exp( -10*ti ) + 3*exp( -4*ti );
    end
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

   x    = varargin{2};
   y    = varargin{3};
   fiel = 0;
   giel = zeros( 6, 1 );
   Hiel = zeros( 6, 6 );
   for iel = 1:13
      ti   = iel / 10;
      riel = x(3) * exp(-x(1)*ti) - x(4)*exp(-x(2)*ti) + x(6)*exp(-x(5)*ti) - y(iel);
      fiel = fiel + riel^2;
      if ( nargout > 1 )
         Jiel = [ -ti*x(3)*exp(-x(1)*ti);  ti*x(4)*exp(-x(2)*ti);  exp(-x(1)*ti);
	            -exp(-x(2)*ti)      ; -ti*x(6)*exp(-x(5)*ti);  exp(-x(5)*ti) ];
	 giel = giel + 2 * Jiel * riel;
         if ( nargout > 2 )
	    Hi =   [ ti^2*x(3)*exp(-x(1)*ti),    0, -ti*exp(-x(1)*ti), 0, 0, 0;
	             0,   -ti^2*x(4)*exp(-x(2)*ti), 0,  ti*exp(-x(2)*ti), 0, 0;
		    -ti*exp(-x(1)*ti),             0,           0,     0, 0, 0;
		     0,               ti*exp(-x(2)*ti),         0,     0, 0, 0;
		     0, 0, 0, 0, ti^2*x(6)*exp(-x(5)*ti), -ti * exp(-x(5)*ti);
		     0, 0, 0, 0,          -ti*exp(-x(5)*ti),    0              ];
            Hiel = Hiel + 2 * ( Jiel * Jiel.' + riel * Hi );
	 end
      end
   end
   varargout{1}= fiel;
   if ( nargout > 1 )
      varargout{2} = giel;
      if ( nargout > 2 )
         varargout{3} = Hiel;
      end
   end

end

return

end