%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = arwhead( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   A variable-dimension quartic problem whose Hessian is an
%   arrow-head (downwards) with diagonal central part and border-width of 1.
%
%   Source: Problem 55 in
%      A.R. Conn, N.I.M. Gould, M. Lescrenier and Ph.L. Toint,
%      "Performance of a multifrontal scheme for partially separable
%      optimization",
%      Report 88/4, Dept of Mathematics, FUNDP (Namur, B), 1988.
%
%   If the dimension is unspecified, the default n = 10 is chosen.
%
%   Ph. Toint 22 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'arwhead';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 1 )
         disp( [ ' ERROR in arwhead: n = ', int2str(n),' but should satisfy n >= 1!' ] )
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
   varargout{8} = 'OUR2-AN-V-0';                % class

case 'cpsstr'

   n     = varargin{1};
   eldom = cell( n-1, 1 );
   for iel = 1:n-1
      eldom{ iel } = [ iel n ];
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

   x    = varargin{2};
   varargout{1} = 3 - 4 * x(1) + ( x(1)^2 + x(2)^2 )^2;
   if ( nargout > 1 )
      varargout{2} = [ -4 + 4 * ( x(1)^2 + x(2)^2 )*x(1);
	                4 * ( x(1)^2 + x(2)^2 )*x(2) ];
      if ( nargout > 2 )
	 varargout{3} = [  12 * x(1)^2 + 4 * x(2)^2,  8 * x(1)*x(2)           ;
	                    8 * x(1)*x(2)          ,  4 * x(1)^2 + 12 * x(2)^2 ];
      end
   end
   
end

return

end