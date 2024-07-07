%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = ncb20c( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   A banded problem with semi-bandwidth 20.  This problem exhibits frequent
%   negative curvature in the exact Hessian. This is a simplified version of
%   NCB20 where the paramter COND has been reduced to 1e2 (from 1e4).
%
%   Source:
%   Ph. Toint, private communication, 1992.
%
%   OPM input: Ph. Toint, 27 VII 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'ncb20c';
problem = str2func( pname );
p       = 20;
ny      = 10;
cond    = 1.0e2;         % instead of 1.0e4

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < p+ny )
         disp( [ ' ERROR in ncb20c: n = ', int2str(n), ' is less than ', num2str(p+ny),'!' ] )
	 return
      end
   else
      n = p+ny;
   end
   varargout{1} = [ zeros(n-ny,1); ones(ny,1)];   % x0
   varargout{2} = '???';                          % fstar
   varargout{3} = '';                             % xtype
   varargout{4} = [];                             % xlower
   varargout{5} = [];                             % xupper
   varargout{6} = [];                             % clower
   varargout{7} = [];                             % cupper
   varargout{8} = 'OUR2-AN-V-0';                  % class

case 'cpsstr'

   n = varargin{1};
   eldom = cell( n+1, 1 );
   for iel = 1:n+1
      if ( iel <= n-p )
         eldom{iel} = [ iel:iel+p-1 ];               % group O involving element E and S
      elseif ( iel <= n )
         eldom{iel} = [ iel ];                       % group O involving element S
      else
         eldom{iel} = [ [ 1:2*ny] [n-ny+1:n] ]; % group W involving elements L
      end
   end
   cpsstr.name  = pname;
   cpsstr.eldom = eldom;
   cpsstr.param = { n, p, ny };
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
   niel  = length( x );
   n     = varargin{3};
   p     = varargin{4};
   ny    = varargin{5};
   riel = 2;
   switch ( nargout )
   case 1
      if ( iel <= n-p )
         riel = riel - (4/p)*sum(x);
         fe   = ncb20c( 'bp', x );
         fs   = ncb20c( 'qr', x(1) );
         riel = riel + 10*fe/iel + fs;
      elseif ( iel <= n )
         fs   = ncb20c( 'qr', x(1) );
         riel = riel + fs; 
      else
         for j = 1:ny
            flj  = ncb20c( '3p', [ x(j) x(ny+j) x(ny+ny+j) ] );
            riel = riel + cond * flj;
	 end
      end
      varargout{1} = riel;
   case 2
      Jiel = zeros( niel, 1 );
      if ( iel <= n-p )
         riel = riel - (4/p)*sum(x);
         [ fe, ge ]  = ncb20c( 'bp', x );
         [ fs, gs ]  = ncb20c( 'qr', x(1) );
         riel = riel + 10*fe/iel + fs;
	 Jiel    = Jiel - (4/p)*ones(niel,1) + 10*ge/iel;
	 Jiel(1) = Jiel(1) + gs;
      elseif ( iel <= n )
         [ fs, gs ]  = ncb20c( 'qr', x(1) );
         riel = riel + fs;
	 Jiel = gs;
      else
         for j = 1:ny
	    ii = [ j ny+j ny+ny+j ];
            [ flj, glj ] = ncb20c( '3p', x(ii) );
            riel     = riel     + cond * flj;
	    Jiel(ii) = Jiel(ii) + cond * glj;
	 end
      end
      varargout{1} = riel;
      varargout{2} = Jiel;
   case 3
      Jiel = zeros( niel, 1 );
      Hiel = zeros( niel, niel );
      if ( iel <= n-p )
         riel = riel- (4/p)*sum(x);
         [ fe, ge, He ]  = ncb20c( 'bp', x );
         [ fs, gs, Hs ]  = ncb20c( 'qr', x(1) );
         riel      = riel + 10*fe/iel + fs;
	 Jiel      = Jiel - (4/p)*ones(niel,1) + 10*ge/iel;
	 Jiel(1)   = Jiel(1)   + gs;
         Hiel      = Hiel      + 10*He/iel;
	 Hiel(1,1) = Hiel(1,1) + Hs;
      elseif ( iel <= n )
         [ fs, gs, Hs ]  = ncb20c( 'qr', x(1) );
         riel = riel + fs;
	 Jiel = gs;
	 Hiel = Hs;
      else
         for j = 1:ny
	    ii = [ j ny+j ny+ny+j ];
            [ flj, glj, Hlj ] = ncb20c( '3p', x(ii) );
            riel        = riel        + cond * flj;
	    Jiel(ii)    = Jiel(ii)    + cond * glj;
	    Hiel(ii,ii) = Hiel(ii,ii) + cond * Hlj;
	 end
      end
      varargout{1} = riel;
      varargout{2} = Jiel;
      varargout{3} = Hiel;
   end

%[iel norm(Jiel) ] %D

case 'bp'

   x   = varargin{1};
   n   = length( x );
   d   = zeros( n, 1 );
   y   = zeros( n, 1 );
   dy  = zeros( n, 1 );
   d2y = zeros( n, 1 );
   for i = 1:n
      d(i)   = 1 + x(i)^2;
      y(i)   = x(i) / d(i);
      dy(i)  = (1-2*x(i)^2/d(i))/d(i);
      d2y(i) = (8*x(i)^3/d(i)-6*x(i))/ d(i)^2;
   end
   sumy = sum(y);
   varargout{1} = sumy^2;
   if ( nargout > 1 )
      varargout{2} = 2*sumy*dy;
      if ( nargout > 2 )
         varargout{3} = 2*(dy*dy'+sumy*diag(d2y));
      end
   end

case 'qr'

   x = varargin{1};
   varargout{1} = x^4;
   if ( nargout > 1 )
      varargout{2} = 4*x^3;
      if ( nargout > 2 )
         varargout{3} = 12*x^2;
      end
   end

case '3p'

   x = varargin{1};
   varargout{1} = x(1)*x(2)*x(3) + 2*x(3)^2;
   if ( nargout > 1 )
      varargout{2} = [ x(2)*x(3); x(1)*x(3); x(1)*x(2)+4*x(3) ];
      if ( nargout > 2 )
         varargout{3} = [  0    x(3)   x(2);
	                  x(3)   0     x(1);
			  x(2)  x(1)    4   ];
      end
   end
   
end

return

end

