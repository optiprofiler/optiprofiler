%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = powr( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The Power problem by Oren.
%
%   Source:
%      S.S. Oren,
%      Self-scaling variable metric algorithms,
%      Part II: implementation and experiments"
%      Management Science 20(5):863-874, 1974.
%   Also problem 179 (p. 83) in
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   The problem name i shortened to "powr" in order to avoid conflict
%   with the built-in MATLAB function.
%
%   If the dimension is unspecified, the default n = 10 is chosen.
%
%   Ph. Toint 22 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'powr';

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 1 )
         disp( [ ' ERROR in powr: n = ', int2str(n), ' but should be > 0!' ] )
      end
   else
      n = 10;
   end
   varargout{1} = ones( n, 1 );                 % x0
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AN-V-0';                % class

case 'cpsstr'

   cpsstr.name  = pname;
   cpsstr.eldom = { [ 1:varargin{1} ] }; 
cpsstr.param = {};
   varargout{1} = cpsstr;

case 'objf'   % varargout = [ f, g, H ]

   x = varargin{1};
   varargout = opm_eval_cpsf( pname, 'elobjf', x, {[1:length(x)]}, nargout, {} );

case 'elobjf'   % varargout = [ f, g, H ]

   x   = varargin{2};
   n   = length(x);
   r   = [ 1:n ]*(x.^2);
   varargout{1} = r^2;
   if ( nargout > 1 )
      J = 2*[1:n].'.*x;
      varargout{2} = 2*J*r;
      if ( nargout > 2 )
         H = 2*diag([1:n]);
         varargout{3} = 2*(J*J.'+r*H);
      end
   end

end

return

end