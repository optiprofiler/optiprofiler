%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = penalty2( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   A penalty function arising from
%   min{ sum_1^n(x_i-1)^2 subject to sum_1^n x_i^2 = 1/4}.
%   At the solution, the dense Hessian has n-2 eigenvalues of order a
%   and two of order 1. The values of a = 1e-5; b = 1 are given by
%   Moré et al., while Buckley proposed different values.
%
%   Source: problem 24 in
%      J.J. More, B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%   Also problem 110 (p. 80) in
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   If the dimension is unspecified, the default n = 10 is chosen.
%
%   Ph. Toint 22 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'penalty2';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 1 )
         disp( [ ' ERROR in penalty2: n = ', int2str(n), ' must be > 1!' ] )
      end
   else
      n = 10;
   end
   varargout{1} = 0.5*ones( n, 1 );             % x0
   varargout{2} = 'unknown';                    % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AN-V-0';                % class

case 'cpsstr'

   n     = varargin{1};
   eldom = cell( 2*n, 1 );
   ie    = 1;
   for iel = 1:n-1
      eldom{ ie }   = [ iel iel+1 ];
      eldom{ ie+1 } = [ iel ];
      ie = ie + 2;
   end
   eldom{ 2*n-1 } = [ 1:n ];
   eldom{ 2*n   } = [ 1 ];
   cpsstr.name  = pname;
   cpsstr.eldom = eldom;
   cpsstr.param = { n };
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

   iel  = varargin{1};
   x    = varargin{2};
   n    = varargin{3};
   a    = 1e-5; b = 1;   % The Moré et al. values
   if ( iel <= 2*n-2 )
      switch ( mod( iel, 2 ) )
      case 1
         i    = (iel+1)/2;
	 ci   = exp(i/10) + exp((i-1)/10);
	 eim1 = exp(x(1)/10);
	 ei   = exp(x(2)/10);
	 riel = eim1+ei-ci;
         varargout{1} = a*riel^2;
         if ( nargout > 1 )
            Jiel         = 0.1 * [ eim1; ei ];
	    varargout{2} = 2*a*Jiel*riel;
            if ( nargout > 2 )
	       Hiel = 0.01 * [ eim1, 0;
	                        0  , ei ];
	       varargout{3} = 2*a*(Jiel*Jiel.'+riel*Hiel);
	    end
	 end
      case 0
%         i    = iel / 2 + 1;
	 ci   = exp(-1/10);
	 ei   = exp(x(1)/10);
	 riel = ei-ci;
         varargout{1} = a*riel^2;
         if ( nargout > 1 )
	    Jiel = 0.1 * ei;
	    varargout{2} = 2*a*Jiel*riel;
            if ( nargout > 2 )
	       Hiel = 0.01*ei;
	       varargout{3} = 2*a*(Jiel*Jiel.'+riel*Hiel);
	    end
	 end
      end
   elseif( iel == 2*n-1 )
      e1n  = [ 1:n ]';
      w    = e1n(n:-1:1);
      riel = b*(w'*(x.^2)-1);
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = 2*w.*x;
	 varargout{2} = 2*b*Jiel*riel;
         if ( nargout > 2 )
	    Hiel = 2*diag(w);
	    varargout{3} = 2*b*(Jiel*Jiel.'+riel*Hiel);
	 end
      end
   else
      riel = x(1) - 1/5;
      varargout{1} = b*riel^2;
      if ( nargout > 1 )
         varargout{2} = 2*b*riel;
         if ( nargout > 2 )
	    varargout{3} = 2*b;
         end
      end
   end
end

return

end