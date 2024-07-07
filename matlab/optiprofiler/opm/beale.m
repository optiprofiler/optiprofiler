%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = beale( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Beale problem in 2 variables.
%   This function is a nonlinear least squares with 15 groups.  
%
%   Source: Problem 5 in
%      J.J. More', B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%   Also problem 89 in
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 17 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   n            = 2;
   varargout{1} = ones( n, 1 );                 % x0
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AN-2-0';                % class

case 'cpsstr'

  cpsstr.name  = 'beale';
  cpsstr.eldom = { [ 1 2 ] };
  cpsstr.param = {};
  varargout{1} = cpsstr;
  
case 'objf'   % varargout = [ f, g, H ]

   x         = varargin{1};
   varargout = opm_eval_cpsf( 'beale', 'elobjf', x, { [ 1 2 ] }, nargout, {} );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

%  iel  = varargin{1};
   x    = varargin{2};
   r1 = 1.5   - x(1) * ( 1 - x(2)   );
   r2 = 2.25  - x(1) * ( 1 - x(2)^2 );
   r3 = 2.625 - x(1) * ( 1 - x(2)^3 );
   varargout{1} = r1^2 + r2^2 + r3^2;
   if ( nargout > 1 )
      J1 = [ x(2)-1   ; x(1) ];
      J2 = [ x(2)^2-1 ; 2 * x(1) * x(2) ];
      J3 = [ x(2)^3-1 ; 3 * x(1) * x(2)^2 ];
      varargout{2} = 2 * ( J1 * r1 + J2 * r2 + J3 * r3 );
      if ( nargout > 2 )
         H1 = [ 0, 1 ;
	        1, 0 ];
	 H2 = [   0   , 2*x(2);
	        2*x(2), 2*x(1)];
	 H3 = [    0    , 3*x(2)^2;
	        3*x(2)^2, 6*x(1)*x(2) ];
         varargout{3} = 2 * ( J1*J1.' + r1*H1 + J2*J2.' + r2*H2 +J3*J3.' + r3*H3 );
      end
   end
end

return

end