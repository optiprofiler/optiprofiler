%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = cosine( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Another function with nontrivial groups and
%   repetitious elements.
%   NB: scaled version of COSINE.  The original version uses scal = 12.
%
%   Source:
%      N. Gould, private communication, 1997
%
%   If the dimension is unspecified, the default n = 10 is chosen.
%
%   Ph. Toint and S. Gratton, 24 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'cosine';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 1 )
         disp( [ ' ERROR in cosine: n = ', int2str(n),' should be > 1!' ] )
      end
   else
      n = 10;
   end
   scal  = 1;                                   % the unscaled version
   varargout{1} = ones( n, 1);                  % x0
   for i = 1:n
      pip  = exp( scal* i /(n-1) );
      varargout{1}(i) = 1 / pip;
   end
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'QUR2-AY-V-0';                % class

case 'cpsstr'

   n   = varargin{1};
   for iel = 1:n-1
      eldom{ iel } = [ iel iel+1 ];
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

   x    = varargin{2};
   fiel = cos( x(1)^2 - 0.5*x(2) );
   varargout{1} = fiel;
   if ( nargout > 1 )
      aux  = sin( x(1)^2 - 0.5*x(2) );
      varargout{2} = - aux*[2*x(1);-0.5];
      if ( nargout > 2 )
         varargout{3} = [ -fiel*(2*x(1))^2-2*aux , x(1)*fiel;
                          x(1)*fiel       ,  -0.25*fiel          ];
      end
   end

return

end