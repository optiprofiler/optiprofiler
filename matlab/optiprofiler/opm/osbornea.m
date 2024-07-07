%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = osbornea( action, varargin )

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

pname   = 'osbornea';
problem = str2func( pname );
switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 0.5; 1.5; -1; 0.01; 0.02 ]; % x0
   varargout{2} = 5.46489e-5;                   % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AN-5-0';                % class

case 'cpsstr'

   cpsstr.name  = pname;
   cpsstr.eldom = { [ 1:5 ] };
   y  = [ 0.844 0.908 0.932 0.936 0.925 0.908 0.881 0.850 0.818 0.784 0.751 0.718 0.685 ...
          0.658 0.628 0.603 0.580 0.558 0.538 0.522 0.506 0.490 0.478 0.467 0.457 0.448 ...
	  0.438 0.431 0.424 0.420 0.414 0.411 0.406 ];
   cpsstr.param = { y };
   varargout{1} = cpsstr;
   
case 'objf'   % varargout = [ f, g, H ]

   x  = varargin{1};
   if ( nargin > 2 && isstruct( varargin{2} ) && isfield( varargin{2}, 'eldom' ) );
      cpsstr = varargin{2};
   else
      cpsstr = problem( 'cpsstr', length(x) );
   end
   varargout = opm_eval_cpsf( pname, 'elobjf', x, cpsstr.eldom, nargout, cpsstr.param{:} );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   x  = varargin{2};
   y  = varargin{3};
   varargout{1} = 0;
   if ( nargout > 1 )
      varargout{2} = zeros( 5, 1 );
      if ( nargout > 2 )
         varargout{3} = zeros( 5, 5 );
      end
   end
   for i = 1:33
      ti = 10*(i-1);
      e4 = exp(-ti*x(4));
      e5 = exp(-ti*x(5));
      ri = x(1) + x(2)*e4 + x(3)*e5 - y(i);
      varargout{1} = varargout{1} + ri^2;
      if ( nargout > 1 )
         Ji = [ 1; e4; e5; -ti*x(2)*e4; -ti*x(3)*e5 ];
	 varargout{2}= varargout{2} + 2*Ji*ri;
         if ( nargout > 2 ) 
	    Hi = [ 0,   0   ,   0    ,      0      ,      0     ;
	           0,   0   ,   0    ,    -ti*e4   ,      0     ;
		   0,   0   ,   0    ,      0      ,   -ti*e5   ;
		   0, -ti*e4,   0    , x(2)*ti^2*e4,      0     ;
		   0,   0   , -ti*e5 ,      0      , x(3)*ti^2*e5 ];
            varargout{3} = varargout{3} +2*(Ji*Ji.' + ri*Hi );
	 end
      end
   end
end

return

end