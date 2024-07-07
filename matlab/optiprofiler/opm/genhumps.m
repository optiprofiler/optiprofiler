%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = genhumps( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   A multi-dimensional variant of HUMPS, a two dimensional function 
%   with a lot of humps. 
%
%   The problem is nonconvex.
%
%   Source: 
%   Ph. Toint, private communication, 1997.
%
%   If the dimension is unspecified, the default n = 12 is chosen.
%
%   Ph. Toint and S. Gratton,  22 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'genhumps';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 2 )
         disp( [ ' ERROR in genhumps: n = ', int2str(n),' should be > 1!' ] )
      end
   else
      n = 10;
   end
   varargout{1} = [ -506; -506.2*ones( n-1, 1) ]; % x0
   varargout{2} = 0;                              % fstar
   varargout{3} = '';                             % xtype
   varargout{4} = [];                             % xlower
   varargout{5} = [];                             % xupper
   varargout{6} = [];                             % clower
   varargout{7} = [];                             % cupper
   varargout{8} = 'SUR2-AY-V-0';                  % class

case 'cpsstr'

   n   = varargin{1};
   for iel = 1:n-1
      eldom{ iel } = [ iel iel+1 ];
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

   x = varargin{2};
   varargout{1} = sin(20*x(1))^2*sin(20*x(2))^2+1/20*(x(1)^2+x(2)^2) ;
   if ( nargout > 1 )
      varargout{2} = [ 20*sin(40*x(1))*sin(20*x(2))^2+1/10*x(1);
                       20*sin(40*x(2))*sin(20*x(1))^2+1/10*x(2) ];
      if ( nargout > 2 )
         varargout{3} = [ 800*cos(40*x(1))*sin(20*x(2))^2+1/10, 400*sin(40*x(1))*sin(40*x(2));
                          400*sin(40*x(1))*sin(40*x(2))       , 800*cos(40*x(2))*sin(20*x(1))^2+1/10];
      end
   end

return

end