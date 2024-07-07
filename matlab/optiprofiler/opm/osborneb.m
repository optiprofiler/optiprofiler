%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = osborneb( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Osborne first problem in 5 variables.
%   This function  is a nonlinear least squares with 33 groups.
%
%   Source:  Problem 17 in
%      J.J. More', B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%   Also problem 32 (p. 77)
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 24 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'osborneb';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 1.3; 0.65; 0.65; 0.7; 0.6; 3.0; 5.0; 7.0; 2.0; 4.5; 5.5 ]; % x0
   varargout{2} = 4.01377e-2;                   % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AN-11-0';               % class

case 'cpsstr'

   cpsstr.name  = pname;
   cpsstr.eldom = { [ 1:11 ] };
   cpsstr.param = { [ 1.366 1.191 1.112 1.013 0.991 0.885 0.831 0.847 0.786 0.725 ...
                      0.746 0.679 0.608 0.655 0.616 0.606 0.602 0.626 0.651 0.724 ...
 	              0.649 0.649 0.694 0.644 0.624 0.661 0.612 0.558 0.533 0.495 ...
	              0.500 0.423 0.395 0.375 0.372 0.391 0.396 0.405 0.428 0.429 ...
	              0.523 0.562 0.607 0.653 0.672 0.708 0.633 0.668 0.645 0.632 ...
	              0.591 0.559 0.597 0.625 0.739 0.710 0.729 0.720 0.636 0.581 ...
	              0.428 0.292 0.162 0.098 0.054 ] };
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

   x  = varargin{2};
   y  = varargin{3};
   varargout{1} = 0;
   if ( nargout > 1 )
      varargout{2} = zeros( 11, 1 );
      if ( nargout > 2 )
         varargout{3} = zeros( 11, 11 );
      end
   end
   for i = 1:65
      ti   = (i-1)/10;
      t9   = ti-x(9);
      t10  = ti-x(10);
      t11  = ti-x(11);
      e5   = exp(-ti*x(5));
      e69  = exp( -x(6)*t9^2);
      e710 = exp( -x(7)*t10^2);
      e811 = exp( -x(8)*t11^2);
      ri = x(1)*e5 + x(2)*e69 + x(3)*e710 +x(4)*e811 - y(i);
      varargout{1} = varargout{1} + ri^2;
      if ( nargout > 1 )
         Ji = [ e5; e69; e710; e811; -ti*x(5)*e5;
	        -x(2)*t9^2*e69; -x(3)*t10^2*e710; -x(4)*t11^2*e811;
		2*x(2)*x(6)*t9*e69; 2*x(3)*x(7)*t10*e710; 2*x(4)*x(8)*t11*e811 ];
	 varargout{2}= varargout{2} + 2*Ji*ri;
         if ( nargout > 2 )
	    Hi        = zeros( 11, 11 );
	    Hi(1,5)   = -ti*e5;
	    Hi(5,1)   = Hi(1,5);
	    Hi(2,6)   = -(ti- x(9))^2*e69;
	    Hi(6,2)   = Hi(2,6);
	    Hi(2,9)   = 2*x(6)*t9*e69;
	    Hi(9,2)   = Hi(2,9);
	    Hi(3,7)   = -(ti- x(10))^2*e710;
	    Hi(7,3)   = Hi(3,7);
	    Hi(3,10)  = 2*x(7)*t10*e710;
	    Hi(10,3)  = Hi(3,10);
	    Hi(4,8)   = -(ti- x(11))^2*e811;
	    Hi(8,4)   = Hi(4,8);
	    Hi(4,11)  =  2*x(8)*t11*e811;
	    Hi(11,4)  = Hi(4,11);
	    Hi(5,5)   = e5*(-ti*+ti^2);
	    Hi(6,6)   = e69*x(2)*t9^4;
	    Hi(6,9)   = 2*x(2)*e69*(t9-x(6)*t9^3);
	    Hi(9,6)   = Hi(6,9);
	    Hi(7,7)   = x(3)*t10^4*e710;
	    Hi(7,10)  = 2*x(3)*e710*(t10-x(7)*t10^3);
	    Hi(10,7)  = Hi(7,10);
	    Hi(8,8)   = x(4)*t11^4*e811;
	    Hi(8,11)  = 2*x(4)*e811*(t11-x(8)*t11^3);
	    Hi(11,8)  = Hi(8,11);
	    Hi(9,9)   = 2*x(2)*x(6)*e69 *(-1+2*x(6)*t9^2);
	    Hi(10,10) = 2*x(3)*x(7)*e710*(-1+2*x(7)*t10^2);
	    Hi(11,11) = 2*x(4)*x(8)*e811*(-1+2*x(8)*t11^2);
            varargout{3} = varargout{3} +2*(Ji*Ji.' + ri*Hi );
	 end
      end
   end
end

return

end