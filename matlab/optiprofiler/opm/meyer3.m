%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = meyer3( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Meyer's nonlinear least-squares problem in 3 variables.
%   A problem arising in the analysis of the resistance of a thermistor.
%
%   Source:  Problem 28 in
%      J.J. More', B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%   Also problem 29 (p. 73)
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 24 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'meyer3';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 0.02; 4000; 250 ];          % x0
   varargout{2} = 87.9458;                      % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-RN-3-0';                % class

case 'cpsstr'

   cpsstr.name  = pname;
   cpsstr.eldom = { [ 1 2 3 ] }; 
   cpsstr.param = { [ 34780 28610 23650 19630 16370 13720 11540 9744 8261 7030 6005 5147 4427 3820 3307 2872 ] };
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
      varargout{2} = zeros( 3, 1 );
      if ( nargout > 2 )
         varargout{3} = zeros( 3, 3 );
      end
   end
   for i = 1:16
      ti    = 45 + 5 * i;
      tipx3 = ti + x(3);
      a     = x(2)/tipx3;
      ea    = exp( a );
      ri    = x(1) * ea - y(i);
      varargout{1} = varargout{1} + ri^2;
      if ( nargout > 1 )
         Ji = ea * [ 1 ; x(1)/tipx3; -x(1)*x(2)/tipx3^2 ];
	 varargout{2}= varargout{2} + 2 * Ji * ri;
         if ( nargout > 2 )
	    Hi        = ea * [    0   ,       1/tipx3      ,    -x(2)/tipx3^2    ;
                               1/tipx3,    x(1)/tipx3^2    ,  -x(1)*(a+1)/tipx3^2;
	                      -a/tipx3, -x(1)*(a+1)/tipx3^2, x(1)*a*(a+2)/tipx3^2];
            varargout{3} = varargout{3} + 2 * ( Ji * Ji.' + ri * Hi );
	 end
      end
   end
end

return

end