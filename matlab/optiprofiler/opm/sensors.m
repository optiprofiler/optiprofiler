%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = sensors( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   A problem arising from two-dimensional optimal sensor placement
%
%   Source:
%   H. Zhang and X. Wang,
%   "Optimal sensor placement",
%   SIAM Review, vol. 35, p. 641, 1993.
%
%   SIF input: Nick Gould, June 1994,
%              OPM with a correction: Ph. Toint 26 VII 2021
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'sensors';
problem = str2func( pname );

switch( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 2 )
         disp( [ ' ERROR in sensors: n = ', int2str(n), ' is less than 2!' ] )
      end
   else
      n = 10;
   end
   varargout{1} = [1:n]'/n;                       % x0
   varargout{2} = '???';                          % fstar
   varargout{3} = '';                             % xtype
   varargout{4} = [];                             % xlower
   varargout{5} = [];                             % xupper
   varargout{6} = [];                             % clower
   varargout{7} = [];                             % cupper
   varargout{8} = 'OUR2-AN-V-0';                  % class

case 'cpsstr'

   n   = varargin{1};
   nel = n;
   eldom = cell( nel, 1 );
   for iel = 1:n
      eldom{iel} = [1:n];
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

   iel  = varargin{1};
   x    = varargin{2};
   n    = length(x);
   riel = 0;
   switch( nargout )
   case 1
      for i = 1:n
         fe  = sensors( 'sinfun', [ x(i); x(iel) ] );
         riel = riel + fe;
      end
      varargout{1} = -riel^2;
   case 2
      Jiel = zeros( n, 1 );
      for i = 1:n
         ii = [ i iel ];
         [ fe, ge ] = sensors( 'sinfun', x(ii) );
         riel       = riel     + fe;
	 Jiel(ii)   = Jiel(ii) + ge;
      end
      varargout{1} = -riel^2;
      varargout{2} = -2*Jiel*riel;
   case 3
      Jiel = zeros( n, 1 );
      Hiel = zeros( n, n );
      for i = 1:n
         ii = [ i iel ];
         [ fe, ge, He ] = sensors( 'sinfun', x(ii) );
         riel           = riel        + fe;
	 Jiel(ii)       = Jiel(ii)    + ge;
	 Hiel(ii,ii)    = Hiel(ii,ii) + He;
      end
      varargout{1} = -riel^2;
      varargout{2} = -2*Jiel*riel;
      varargout{3} = -2*(Jiel*Jiel' + riel*Hiel);
   end

case 'sinfun'

   thetai  = varargin{1}(1);
   thetaj  = varargin{1}(2);
   timj    = thetai - thetaj;
   if ( timj ~= 0 )%D
      si      = sin( thetai );
      sj      = sin( thetaj );
      simj    = sin( timj );
      varargout{1} = si*sj*simj;
      if ( nargout > 1 )
         ci     = cos( thetai );
         cj     = cos( thetaj );
         cimj   = cos( timj );
         cjsimj = cj*simj - sj*cimj;
         varargout{2} = [ sj*(ci*simj+si*cimj); si*cjsimj ];
         if ( nargout > 2 )
            varargout{3} = [ 2*sj*(ci*cimj-si*simj)          ci*cjsimj+si*(cj*cimj+sj*simj)
                             ci*cjsimj+si*(cj*cimj+sj*simj)  -si*(sj*simj+cj*cimj+cj*cimj+sj*simj) ];
         end
      end
   else
      varargout{1} = 0;
      varargout{2} = zeros( 2, 1 );
      varargout{3} = zeros( 2, 2 );
   end

end

return

end
