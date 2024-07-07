%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = watson( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Watson problem in varaible dimension ( 2 <= n <= 31 ).
%   This function is a nonlinear least squares with 31 groups. 
%
%   Source:  problem 20 in
%      J.J. More', B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%   Also problem 128 (p. 100) in 
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   If the dimension is unspecified, the default n = 12 is chosen
%   (as in CUTest).
%
%   Ph. Toint 22 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'watson';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 2 || n > 31 )
         disp( [ ' ERROR in watson: n = ', int2str(n),' but should satisfy 2 <= n <= 31!' ] )
      end
   else
      n = 12;
   end
   varargout{1} = zeros( n, 1 );                % x0
   varargout{2} = 'unknown';                    % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AN-V-0';                % class

case 'cpsstr'

   cpsstr.name  = pname;
   cpsstr.eldom = { [1:varargin{1}] };
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

   x = varargin{2};
   n = length( x );
   varargout{1} = 0;
   if ( nargout > 1 )
      varargout{2} = zeros( n, 1 );
      if ( nargout > 2 )
         varargout{3} = sparse( n, n );
      end
   end
   for i=1:29
      w1 = zeros( n, 1 );
      w2 = zeros( n, 1 );
      for j = 1:n
         w1(j) = (j-1)*(i/29)^(j-2);
	 w2(j) = (i/29)^(j-1);
      end
      riel = w1'*x - (w2'*x)^2 - 1;
      varargout{1} = varargout{1} + riel^2;
      if ( nargout > 1 )
         Jiel = w1-2*(w2'*x)*w2;
	 varargout{2} = varargout{2}+2*Jiel*riel;
         if ( nargout > 2 )
	    Hiel = -2*w2*w2.';
	    varargout{3} = varargout{3} + 2*(Jiel*Jiel.'+riel*Hiel);
	 end
      end
   end
   varargout{1} = varargout{1}+x(1)^2;
   if ( nargout > 1 )
      varargout{2} = varargout{2} + [ 2*x(1); zeros( n-1, 1 ) ];
      if ( nargout > 2 )
         Hiel         = sparse(n,n);
	 Hiel( 1, 1 ) = 2;
         varargout{3} = varargout{3} + Hiel;
      end
   end
   riel = x(2)-x(1)^2-1;
   varargout{1} = varargout{1}+riel^2;
   if ( nargout > 1 )
      Jiel = [ -2*x(1); 1; zeros( n-2, 1) ];
      varargout{2} = varargout{2} + 2*Jiel*riel;
      if ( nargout > 2 )
         Hiel         = sparse(n,n);
	 Hiel( 1, 1 ) = -2;
         varargout{3} = varargout{3} + 2*(Jiel*Jiel.'+riel*Hiel);
      end
   end

end

return

end