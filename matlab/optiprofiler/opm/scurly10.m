%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = scurly10( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   A banded function with semi-bandwidth 20 and
%   negative curvature near the starting point.
%   NB: scaled version of CURLY10.
%
%   Source: Nick Gould
%
%   SIF input: Nick Gould, September 1997.
%              OPM: Ph. Toint, 26 VII 2021
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'scurly10';
problem = str2func( pname );
smb     = 10;  % semibandwidth
scl     = 12;  % the range of scales

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < smb )
         disp( [ ' ERROR in scurly10: n = ', int2str(n), ' is less than ', num2str(smb), '!' ] )
      end
   else
      n = 30;
   end
   scale        = exp( [0:n-1]/(n-1)*scl )';
   varargout{1} = 0.0001*scale.*[1:n]'/(n+1);     % x0
   varargout{2} = '???';                          % fstar
   varargout{3} = '';                             % xtype
   varargout{4} = [];                             % xlower
   varargout{5} = [];                             % xupper
   varargout{6} = [];                             % clower
   varargout{7} = [];                             % cupper
   varargout{8} = 'OUR2-AN-V-0';                  % class

case 'cpsstr'

   n     = varargin{1};
   scale = exp( [0:n-1]/(n-1)*scl )';
   eldom = cell( n, 1 );
   for iel = 1:n
      eldom{iel}   = [ iel:min(iel+smb,n) ];
   end
   cpsstr.name  = pname;
   cpsstr.eldom = eldom;
   cpsstr.param = { eldom, scale };
   varargout{1} = cpsstr;

case 'objf'   % varargout = [ f, g, H ]

   x = varargin{1};
   n = length( x );
   if ( nargin > 2 && isstruct( varargin{2} ) && isfield( varargin{2}, 'eldom' ) )
      cpsstr = varargin{2};
   else
      cpsstr = problem( 'cpsstr', n );
   end
   varargout = opm_eval_cpsf( pname, 'elobjf', x, cpsstr.eldom, nargout, cpsstr.param{:} );
   
case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   iel   = varargin{1};
   x     = varargin{2};
   eldom = varargin{3};
   scalx = varargin{4}(eldom{iel});
   x     = x .* scalx;
   switch( nargout )
   case 1
      varargout{1} = scurly10( 'p4', sum(x) );
   case 2
      [ varargout{1}, gsum ] = scurly10( 'p4', sum(x) );
      varargout{2} = gsum * scalx;
   case 3
      [ varargout{1}, gsum, Hsum ] = scurly10( 'p4', sum(x) );
      varargout{2} = gsum * scalx;
      varargout{3} = Hsum * scalx*scalx';
   end
   
case 'p4'

    v  = varargin{1};
    varargout{1} = v * ( v * ( v^2 - 20 ) - 0.1 );
    if ( nargout > 1 )
       varargout{2} = 2 * v * ( 2*v^2 - 20 ) - 0.1;
       if ( nargout > 2 )
          varargout{3} = 12 * v^2 - 2 * 20;
       end
    end

end

return

end
