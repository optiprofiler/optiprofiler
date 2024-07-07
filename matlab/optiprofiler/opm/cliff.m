%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = cliff( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Cliff problem in 2 variables.
%
%   Source: problem 206 in 
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 17 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'cliff';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 0; -1];                     % x0
   varargout{2} = 0.19978661;                   % fstar
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
   ee   = exp(20*(x(1)-x(2)));
   varargout{1} = ((x(1)-3)/100)^2 - (x(1) - x(2)) + ee;
   if ( nargout > 1 )
      varargout{2} = [ 0.0002*(x(1)-3)-1+20*ee; 1-20*ee ];
      if ( nargout > 2 )
	 varargout{3} = [  0.0002+400*ee, -400*ee;
	                   -400*ee       ,  400*ee ];
      end
   end
end

return

end