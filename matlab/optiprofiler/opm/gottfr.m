%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = gottfr( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   A simple two-diemnsional nonlinear least-squares problem.
%
%   Source: problem 208 (p. 56) in
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 20 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'gottfr';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 0.5; 0.5];                  % x0
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AN-2-0';                % class

case 'cpsstr'

   eldom{ 1 } = [ 1 2 ];
   eldom{ 2 } = [ 1 2 ];
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

   iel   = varargin{1};
   x     = varargin{2};
   switch ( iel )
   case 1
      riel = x(1) - 0.1136 * ( x(1) + 3 * x(2) ) * ( 1 - x(1) );
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = [ 1-0.1136*(1-2*x(1)-3*x(2)); -0.1136*3*(1-x(1)) ];
	 varargout{2} = 2 * Jiel * riel;
         if ( nargout > 2 )
	    Hiel = [ 2*0.1136, 3*0.1136;
	             3*0.1136,     0   ];
            varargout{3} = 2 * ( Jiel*Jiel.' + riel*Hiel );
	 end
      end
   case 2
      riel = x(2) + 7.5 * ( 2 * x(1) - x(2) ) * ( 1 - x(2) );
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = [15*(1-x(2)); 1+7.5*(2*x(2)-2*x(1)-1) ];
	 varargout{2} = 2 * Jiel * riel;
         if ( nargout > 2 )
	    Hiel = [   0, -15;
	             -15,  15 ];
            varargout{3} = 2 * ( Jiel*Jiel.' + riel*Hiel );
	 end
      end
   end

return

end