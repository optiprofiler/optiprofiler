%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = box2( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Box problem in 2 variables.
%   This function is a nonlinear least squares with 1O groups.  It is a
%   variation on the box3 problem where x(3) variable is fixed to 1.
%
%   Source: Inspited by problem 12 in
%      J.J. More', B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%   Also problem 11 in 
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 17 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'box2';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 0; 10; 1];                  % x0
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [ -Inf; -Inf; 1 ];            % xlower
   varargout{5} = [  Inf;  Inf; 1 ];            % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'OXR2-AN-2-0';                % class

case 'cpsstr'

   for iel = 1:10
      eldom{ iel } = [ 1 2 3 ];
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

   %   Zero derivatives of fixed variables.

   if ( nargout > 1 )
      varargout{2}(3) = 0;
      if ( nargout > 2 )
         varargout{3}( 3, 1:3 ) = 0;
	 varargout{3}( 1:3, 3 ) = 0;
      end
   end
   
case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   iel  = varargin{1};
   x    = varargin{2};
   ti   = iel / 10;
   riel = exp(-x(1)*ti) - exp(-x(2)*ti) - x(3)*(exp(-ti)-exp(-iel));
   varargout{1} = riel^2;
   if ( nargout > 1 )
      Jiel = [ -ti*exp(-x(1)*ti); ti*exp(-x(2)*ti); -(exp(-ti)-exp(-iel)) ];
      varargout{2} = 2 * Jiel * riel;
      if ( nargout > 2 )
         Hiel = [ ti^2*exp(-x(1)*ti),  0, 0;
	          0, -ti^2*exp(-x(2)*ti), 0;
		  0,         0            0 ];
	 varargout{3} = 2 * ( Jiel*Jiel.' + riel*Hiel );
      end
   end
end

return

end