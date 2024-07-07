%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = sisser( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   A simple unconstrained problem in 2 variables.
%
%   Source:
%      F.S. Sisser,
%      "Elimination of bounds in optimization problems by transforming
%      variables",
%      Mathematical Programming 20:110-121, 1981.
%   Also problem 216 (p. 91) in 
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 22 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'sisser';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 1; 0.1 ];                  % x0
   varargout{2} = 0;                            % fstar
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
   if ( nargin > 2 && isstruct( varargin{2} ) && isfield( varargin{2}, 'eldom' ) )
      cpsstr = varargin{2};
   else
      cpsstr = problem( 'cpsstr', length( x ) );
   end
   varargout = opm_eval_cpsf( pname, 'elobjf', x, cpsstr.eldom, nargout, {} );

case 'elobjf'   % varargout = [ f, g, H ]

   x       = varargin{2};
   varargout{1} = 3*x(1)^4 - 2*x(1)^2*x(2)^2 + 3*x(2)^4;
   if ( nargout > 1 )
      varargout{2} = [ 12*x(1)^3-4*x(1)*x(2)^2; -4*x(1)^2*x(2)+12*x(2)^3 ];
      if ( nargout > 2 )
         varargout{3} = [ 36*x(1)^2-4*x(2)^2, -8*x(1)*x(2);
	                  -8*x(1)*x(2), 36*x(2)^2-4*x(1)^2 ];
      end
   end

end

return

end