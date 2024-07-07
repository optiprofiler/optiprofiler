%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = chandheu( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Chandrasekhar Radiative Transfer H equation, as stated by T. Kelley
%   and expressed in the least-squares sense.
%
%   Source: problem 4 in
%   J.J. More',
%   "A collection of nonlinear model problems"
%   Proceedings of the AMS-SIAM Summer seminar on the Computational
%   Solution of Nonlinear Systems of Equations, Colorado, 1988.
%   Argonne National Laboratory MCS-P60-0289, 1989.
%   SIF input: Ph. Toint, Dec 1989.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'chandheu';
problem = str2func( pname );
switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 10 )
         disp( [ ' ERROR in hydc20ls: n = ', int2str(n), ' is not equal to 57!' ] )
      end
   else
      n = 57;
   end
   varargout{1} = ones(n,1);                      % x0
   varargout{2} = 0;                              % fstar
   varargout{3} = '';                             % xtype
   varargout{4} = [];                             % xlower
   varargout{5} = [];                             % xupper
   varargout{6} = [];                             % clower
   varargout{7} = [];                             % cupper
   varargout{8} = 'SUR2-RN-V-0';                  % class

case 'cpsstr'

   n = varargin{1};        % the number of discretization points

%  Discretization points (X(I)) and weights (W(I)) for the considered
%  quadrature rule on [0,1]

   eldom = cell( n, 1 );
   for iel = 1:n
      eldom{iel} = [ 1:n ];
   end
   
%  The value of the problem parameter C should be in [0,1] for a
%  physically realistic problem, but other values can be used for
%  testing purposes.
%  Unique solution for C=0 and C=1, two for other values.
%  More difficult for C close to 1.

   c = 1;
 
   cpsstr.name  = pname;
   cpsstr.eldom = eldom;
   cpsstr.param = { [1:n]/n, 0.5*c/n };
   varargout{1} = cpsstr;

case 'objf'   % varargout = [ f, g, H ]

   x = varargin{1};
   if ( nargin > 2 && isstruct( varargin{2} ) && isfield( varargin{2}, 'eldom' ) );
      cpsstr = varargin{2};
   else
      cpsstr = problem( 'cpsstr', length(x) );
   end
   varargout = opm_eval_cpsf( pname, 'elobjf', x, cpsstr.eldom, nargout, cpsstr.param{:} );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   iel   = varargin{1};
   x     = varargin{2};
   n     = length( x );
   xx    = varargin{3};
   hcw   = varargin{4};
   riel  = 0;
   Jiel  = zeros( n, 1 );
   Hiel  = sparse( n, n );
   for j = 1:n
      aij  = -xx(iel)*hcw/(xx(iel)+xx(j));
      riel = riel + aij * x(iel) * x(j) + x(iel);
      if ( nargout > 1 )
         if ( iel == j )
	    Jiel(iel) = Jiel(iel) + 2*aij*x(iel) + 1;
	    if ( nargout > 2 )
               Hiel( iel, iel ) = Hiel(iel,iel) + 2*aij;
	    end
	 else
            Jiel(iel) = Jiel(iel) +  aij*x(j) + 1;
	    Jiel(j)   = Jiel(j)   +  aij*x(iel);
	    if ( nargout > 2 )
	       Hiel( iel, j ) = Hiel( iel, j ) + aij;
	       Hiel( j, iel ) = Hiel( j, iel ) + aij;
	    end
	 end
      end
   end
   varargout{1} = riel^2;
   if ( nargout > 1 )
      varargout{2} = 2*Jiel*riel;
      if ( nargout > 2 )
         varargout{3} = 2*( Jiel*Jiel' + riel*Hiel );
      end
   end	 
end

return

end

