%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = powellbs( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Powell's badly scaled problem.
%
%   Source:  Problem 3 in
%      J.J. More', B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%   Also problem 22 (p. 82)
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 24 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'powellbs';

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 0; 1];                      % x0
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

   varargout = opm_eval_cpsf( pname, 'elobjf', varargin{1}, { [ 1 2 ] }, nargout, {} );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   x  = varargin{2};
   a  = 1e4;
   r1 = a * x(1) * x(2) - 1;
   e1 = exp(-x(1));
   e2 = exp(-x(2));
   r2 = e1 + e2 - 1.0001;
   varargout{1} = r1^2 + r2^2;
   if ( nargout > 1 )
      J1 = a * [ x(2); x(1) ];
      J2 = [ -e1; -e2 ];
      varargout{2} = 2 * ( J1 * r1  + J2 * r2 );
      if ( nargout > 2 )
         H1 = a * [ 0, 1;
	            1, 0 ];
         H2 = [ e1, 0;
	        0 , e2 ];
         varargout{3} = 2 * ( J1 * J1.' + r1 * H1 + J2 * J2.' + r2 * H2 );
      end
   end
end

return

end