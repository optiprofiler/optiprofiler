%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = indef( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   A nonconvex problem which has an indefinite Hessian at
%   the starting point.
%
%   Source:
%      Nick Gould, CUTE, Oct 1992.
%
%   If the dimension is unspecified, the default n = 10 is chosen.
%
%   Ph. Toint 22 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'indef';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 2 )
         disp( [ ' ERROR in indef: n = ', int2str(n),' does not satisfy n > 1!' ] )
      end
   else
      n = 10;
   end
   varargout{1} = [ 1:n ]'/(n+1);               % x0
   varargout{2} = 'unknown';                    % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'OUR2-AY-V-0';                % class

case 'cpsstr'

   n       = varargin{1};
   eldom   = cell( n, 1 );
   for iel = 1:n                         %
      if ( iel == 1 || iel == n )
         eldom{ iel } = [ iel ];
      else
         eldom{ iel } = [ 1 iel n ];
      end
   end
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

   iel = varargin{1};
   x   = varargin{2};
   n   = varargin{3};
   scale = 100;
   if ( iel == 1 || iel == n )
      varargout{1} = scale * sin( x(1)/scale ) ;
   else
      aux1 = cos( -x(1)+2*x(2)-x(3) );
      varargout{1} = 0.5 * aux1 + scale * sin( x(2)/scale ) ;
   end
   if ( nargout > 1 )
      if ( iel == 1 || iel == n )
         varargout{2} = cos( x(1)/scale ) ;
      else
         aux2 = sin( -x(1)+2*x(2)-x(3) );
         varargout{2} = [ 0.5*aux2; -aux2 + cos( x(2)/scale ); 0.5*aux2 ];
      end
      if ( nargout > 2 )
         if ( iel == 1 || iel == n )
            varargout{3} = -(1/scale)*sin( x(1)/scale );
         else
            varargout{3} = [ -0.5*aux1,    aux1                             , -0.5*aux1;
                                  aux1, -2*aux1-(1/scale)*sin( x(2)/scale ) ,      aux1;
                             -0.5*aux1,    aux1                             , -0.5*aux1 ];
         end
      end
   end

end

return

end