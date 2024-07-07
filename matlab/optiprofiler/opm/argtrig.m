%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = argtrig( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Variable dimension trigonometric problem
%   This problem is a sum of n least-squares groups, each of
%   which has n+1 nonlinear elements.  Its Hessian matrix is dense.
%
%   Source: problem 26 in
%      J.J. More, B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%   Also problem 101 in
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

pname   = 'argtrig';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 1 )
         disp( [ ' ERROR in argtrig: n = ', int2str(n),' but should satisfy n >= 1!' ] )
      end
   else
      n = 10;
   end
   varargout{1} = ones( n, 1 );                 % x0
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AN-V-0';                % class

case 'cpsstr'

   n     = varargin{1};
   eldom = cell( n, 1 );
   for iel = 1:n
      eldom{ iel } = [ 1:n ];
   end
   cpsstr.name  = pname;
   cpsstr.eldom = eldom;
   cpsstr.param = {};
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
   n    = length( x );
   siel = sin( x(iel) );
   ciel = cos( x(iel) );
   riel  = n - iel * ( 1 - ciel ) - siel;
   for j = 1:n
       riel = riel - cos( x(j) );
   end
   varargout{1} = riel^2;
   if ( nargout > 1 )
      Jiel = zeros( n, 1 );
      Jiel( iel ) = -iel * siel - ciel;
      for j = 1:n
          Jiel( j ) = Jiel( j ) + sin( x(j) );
      end
      varargout{2} = 2 * Jiel * riel;
      if ( nargout > 2 )
         Hiel = zeros( n, n );
         Hiel(iel,iel) = -iel * ciel + siel;
         for j = 1:n
            Hiel(j,j) = Hiel(j,j) + cos( x(j) );
         end
         varargout{3} = 2 * ( Jiel * Jiel.' + riel * Hiel );
      end
   end
   
end

return

end