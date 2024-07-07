%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = nzf1( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The NZF1 problem
%
%   Source:
%      C. Price and Ph. L. Toint,
%      "Exploiting problem structure in pattern-search methods 
%        for unconstrained optimization",
%      Optimization Methods and Software, vol. 21(2), pp. 479-491, 2006.
%
%   The problem dimension must be a multiple of 13.
%   If the dimension is unspecified, the default n = 13 is chosen.
%
%   Ph. Toint 21 VII 2021.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'nzf1';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 13 || round( n/13 ) ~= n/13 )
         disp( [ ' ERROR in nzf1: n = ', int2str(n), ' is not a multiple of 13!' ] )
      end
   else
      n = 13;
   end
   varargout{1} = ones(n,1);                      % x0
   varargout{2} = 0;                              % fstar
   varargout{3} = '';                             % xtype
   varargout{4} = [];                             % xlower
   varargout{5} = [];                             % xupper
   varargout{6} = [];                             % clower
   varargout{7} = [];                             % cupper
   varargout{8} = 'OUR2-AY-V-0';                  % class

case 'cpsstr'

   m     = varargin{1} / 13;
   nel   = 1;
   eldom = cell( m+4, 1 );
   for iel = 1:m
      eldom{ nel }     = [ iel iel+1 iel+2 ];
      eldom{ nel + 1 } = [ iel+1 iel+2 iel+3 iel+4 iel+5 iel+6 ];
      eldom{ nel + 2 } = [ iel+5 iel+7 iel+8 iel+10 ];
      eldom{ nel + 3 } = [ iel+10 iel+11 iel+12 ];
      eldom{ nel + 4 } = [ iel+4 iel+5 iel+9 ];
      nel              = nel + 5;
   end
   for iel = 1:m-1
      eldom{ nel } = [ iel+6 iel+19 ];
      nel          = nel +1;     
   end
   cpsstr.name  = pname;
   cpsstr.eldom = eldom;
   cpsstr.param = { m };
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

   iel = varargin{1};
   x   = varargin{2};
   m   = varargin{3};
   if ( iel <= 5*m )
      switch( mod( iel, 5 ) )
      case 1
         riel = 3*x(1) - 60 + 0.1*(x(2)-x(3))^2;
         varargout{1} = riel^2;
	 if ( nargout > 1 )
	    Jiel = [ 3; 0.2*(x(2)-x(3)); 0.2*(x(3)-x(2)) ];
            varargout{2} = 2 * Jiel * riel;
	    if ( nargout > 2 )
	       Hiel = [ 0     0        0  ;
	                0     0.2     -0.2;
	                0    -0.2      0.2 ];
               varargout{3} = 2 * ( Jiel*Jiel' + riel*Hiel );
	    end
	 end
      case 2
         riel = x(1)^2 + x(2)^2 + x(3)^2*(1+x(3))^2 + x(6) + x(5)/(1+x(4)^2+sin(0.001*x(4)));
         varargout{1} = riel^2;
	 if ( nargout > 1 )
	    Jiel = [ 2*x(1);
	             2*x(2);
		     2*x(3)*(1+x(3))^2+2*x(3)^2*(1+x(3)); 
	            -x(5)*(2*x(4)+0.001*cos(0.001*x(4)))/(1+x(4)^2+sin(0.001*x(4)))^2;
		     1/(1+x(4)^2+sin(0.001*x(4)));
		     1 ];
            varargout{2} = 2 * Jiel * riel;
	    if ( nargout > 2 )
	       Hiel      = zeros(6,6);
	       Hiel(1,1) = 2;
	       Hiel(2,2) = 2;
	       Hiel(3,3) = 2*(1+x(3))*(1+2*x(3))+2*x(3)*(1+2*x(3))+4*x(3)*(1+x(3));
	       Hiel(4,5) = -(2*x(4)+0.001*cos(0.001*x(4)))/(1+x(4)^2+sin(0.001*x(4)))^2;
	       Hiel(5,4) = Hiel(4,5);
               varargout{3} = 2 * ( Jiel*Jiel' + riel*Hiel );
	    end
	 end
      case 3
         riel = x(1) + x(2) - x(3)^2 + x(4);
         varargout{1} = riel^2;
	 if ( nargout > 1 )
	    Jiel = [ 1; 1; -2*x(3); 1 ];
            varargout{2} = 2 * Jiel * riel;
	    if ( nargout > 2 )
	       Hiel         = zeros(4,4);
	       Hiel(3,3)    = -2;
               varargout{3} = 2 * ( Jiel*Jiel' + riel*Hiel );
	    end
	 end
      case 4
         riel = log( 1 + x(1)^2) + x(2) - 5*x(3) + 20;
         varargout{1} = riel^2;
	 if ( nargout > 1 )
	    Jiel = [ 2*x(1)/( 1 + x(1)^2); 1; -5 ];
            varargout{2} = 2 * Jiel * riel;
	    if ( nargout > 2 )
	       Hiel      = zeros(3,3);
	       Hiel(1,1) = 2/( 1 + x(1)^2) - 4*x(1)^2/( 1 + x(1)^2)^2;
               varargout{3} = 2 * ( Jiel*Jiel' + riel*Hiel );
	    end
	 end
      case 0
         riel = x(1) + x(2) + x(2)*x(3) + 10*x(3) - 50;
         varargout{1} = riel^2;
	 if ( nargout > 1 )
	    Jiel = [ 1; 1+x(3); x(2)+10 ];
            varargout{2} = 2 * Jiel * riel;
	    if ( nargout > 2 )
	       Hiel = [  0   0   0;
	                 0   0   1;
			 0   1   0 ];
               varargout{3} = 2 * ( Jiel*Jiel.' + riel*Hiel );
	    end
	 end
      end
   else
      riel = x(1) - x(2);
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = [ 1; -1 ];
         varargout{2} = 2 * Jiel * riel;
         if ( nargout > 2 )
            varargout{3} = 2 * Jiel*Jiel';
         end
      end
   end
end

return

end