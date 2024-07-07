%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = woods( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The extended Woods problem (separable version).
%   This version uses a slightly unorthodox expression of Woods
%   function as a sum of squares (see Buckley)
%
%   Source:  problem 14 in
%      J.J. More', B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%   Also problem 17 (p. 101) in 
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989
%   and problem 27 in 
%      Ph.L. Toint,
%      "Test problems for partially separable optimization and results
%      for the routine PSPMIN",
%      Report 83/4, Department of Mathematics, FUNDP (Namur, B), 1983.
%
%   The problem dimension must be a multiple of 4.
%   If the dimension is unspecified, the default n = 12 is chosen.
%
%   Ph. Toint 22 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'woods';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 4 || round( n/4 ) ~= n/4 )
         disp( [ ' ERROR in woods: n = ', int2str(n), ' is not a multiple of 4!' ] )
      end
   else
      n = 12;
   end
   varargout{1}(1:2:n-1,1) = -3 * ones( n/2, 1 ); % x0
   varargout{1}(2:2:n,1)   = -1 * ones( n/2, 1 ); % x0
   varargout{2} = 0;                              % fstar
   varargout{3} = '';                             % xtype
   varargout{4} = [];                             % xlower
   varargout{5} = [];                             % xupper
   varargout{6} = [];                             % clower
   varargout{7} = [];                             % cupper
   varargout{8} = 'OUR2-AY-V-0';                  % class

case 'cpsstr'

   m  = varargin{1} / 4;
   ie = 1;
   eldom = cell( m+6, 1 );
   for iel = 1:m
      eldom{ ie }     = [ 4*iel-3 4*iel-2 ];
      eldom{ ie + 1 } = [ 4*iel-3 ];
      eldom{ ie + 2 } = [ 4*iel-1 4*iel ];
      eldom{ ie + 3 } = [ 4*iel-1 ];
      eldom{ ie + 4 } = [ 4*iel-2 ];
      eldom{ ie + 5 } = [ 4*iel ];
      eldom{ ie + 6 } = [ 4*iel-2 4*iel ];
      ie              = ie + 7;
   end
   cpsstr.name  = pname;
   cpsstr.eldom = eldom;
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

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   iel = varargin{1};
   x   = varargin{2};
   switch( mod( iel, 7 ) )
   case 1
      riel = 10*( x(2)-x(1)^2 );
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = 10 * [ -2*x(1); 1 ];
	 varargout{2} = 2 * Jiel * riel;
	 if ( nargout > 2 )
	    Hiel = 10 * [ -2, 0;
	                  0, 0 ];
            varargout{3} = 2 * ( Jiel*Jiel.' + riel*Hiel );
	 end
      end
   case { 2, 4 }
      riel = 1- x(1);
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = [ -1 ];
	 varargout{2} = 2 * Jiel * riel;
	 if ( nargout > 2 )
            varargout{3} = 2 * Jiel*Jiel.';
	 end
      end
   case 3
      riel = sqrt(90)*(x(2)-x(1)^2);
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = sqrt(90) * [ -2*x(1); 1 ];
	 varargout{2} = 2 * Jiel * riel;
	 if ( nargout > 2 )
	    Hiel = sqrt(90)*[ -2, 0;
	                       0, 0 ];
            varargout{3} = 2 * ( Jiel*Jiel.' + riel*Hiel );
	 end
      end
   case { 5, 6 }
      riel = sqrt(10.1)*(x(1)-1);
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = sqrt(10.1);
	 varargout{2} = 2 * Jiel * riel;
	 if ( nargout > 2 )
            varargout{3} = 2 * Jiel*Jiel.';
	 end
      end
   case 0
      riel = sqrt(19.8)*(x(1)-1)*(x(2)-1);
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = sqrt(19.8)*[ x(2)-1; x(1)-1 ];
	 varargout{2} = 2 * Jiel * riel;
	 if ( nargout > 2 )
            Hiel = sqrt(19.8)*[ 0, 1;
	                        1, 0 ];
            varargout{3} = 2 * ( Jiel*Jiel.' + riel*Hiel );
	 end
      end
   end
end

return

end