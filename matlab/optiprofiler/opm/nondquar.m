%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = nondquar( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   A nondiagonal quartic test problem.
%
%   This problem has an arrow-head type Hessian with a tridiagonal
%   central part and a border of width 1.
%   The Hessian is singular at the solution.
%
%   Source: problem 57 in
%   A.R. Conn, N.I.M. Gould, M. Lescrenier and Ph.L. Toint,
%   "Performance of a multi-frontal scheme for partially separable
%   optimization"
%   Report 88/4, Dept of Mathematics, FUNDP (Namur, B), 1988.
%
%   n must be even.
%
%   SIF input: Ph. Toint, Dec 1989.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'nondquar';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( mod(n,2) )
         disp( [ ' ERROR in nondquar: n = ', int2str(n), ' is not even!' ] )
      end
   else
      n = 100;
   end
   pones = ones(n/2,1);
   varargout{1}(1:2:n,1) =  pones;
   varargout{1}(2:2:n,1) = -pones;                % x0
   varargout{2} = 0;                              % fstar
   varargout{3} = '';                             % xtype
   varargout{4} = [];                             % xlower
   varargout{5} = [];                             % xupper
   varargout{6} = [];                             % clower
   varargout{7} = [];                             % cupper
   varargout{8} = 'OUR2-AN-V-0';                  % class

case 'cpsstr'

   n     = varargin{1};
   eldom = cell( n, 1 );
   for iel = 1:n-2
      eldom{iel} = [ iel iel+1 n ];
   end
   eldom{n-1}   = [  1  2 ];
   eldom{n}     = [ n-1 n ];
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

   iel   = varargin{1};
   x     = varargin{2};
   n     = varargin{3};
   if ( iel <= n-2 )
      riel = x(1) + x(2) + x(3);
      varargout{1} = riel^4;
      if ( nargout > 1 )
         varargout{2} = 4*riel^3*ones(3,1);
	 if ( nargout > 2 )
	    varargout{3} = 12*riel^2* ones(3,3);
	 end
      end
   else
      riel = x(1) - x(2);
      varargout{1} =  riel^2;
      if ( nargout > 1 )
         varargout{2} = 2* riel * [ 1; -1 ];
	 if ( nargout > 2 )
            varargout{3} = 2*[ 1; -1 ]*[ 1; -1 ]';
	 end
      end
   end

end

return

end
