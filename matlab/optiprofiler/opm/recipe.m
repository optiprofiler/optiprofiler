%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = recipe( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Source: problem 155 (p. 88) in
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 21 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'recipe';

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 2; 5; 1 ];                  % x0
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AN-3-0';                % class

case 'cpsstr'

   cpsstr.name  = pname;
   cpsstr.eldom = { [1:varargin{1}] };
   cpsstr.param = {};
   varargout{1} = cpsstr;

case 'objf'   % varargout = [ f, g, H ]

   x = varargin{1};
   varargout = opm_eval_cpsf( pname, 'elobjf', x, {[1:length(x)]}, nargout, {} );

case { 'elobjf' }   % varargout = [ f, g, H ]

   x     = varargin{2};
   r1    = x(1)-5;
   r2    = x(2);
   x2mx1 = x(2)-x(1);
   r3    = x(3)/x2mx1;
   varargout{1} = r1^2+r2^2+r3^2;
   if ( nargout > 1 )
      J1 = [ 1; 0; 0 ];
      J2 = [ 0; 1; 0 ];
      J3 = [x(3)/x2mx1^2; -x(3)/x2mx1^2; 1/x2mx1 ];
      varargout{2} = 2*(J1*r1+J2*r2+J3*r3);
      if ( nargout > 2 )
         H3 = [ 2*x(3)/x2mx1^3, -2*x(3)/x2mx1^3, 1/x2mx1^2;
	       -2*x(3)/x2mx1^3,  2*x(3)/x2mx1^3,-1/x2mx1^2;
	          1/x2mx1^2   ,    -1/x2mx1^2  ,   0       ];
         varargout{3} = 2*(J1*J1.'+J2*J2.'+J3*J3.'+r3*H3);
      end
   end

end

return

end