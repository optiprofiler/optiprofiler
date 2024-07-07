%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = eigenals( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Solving symmetric eigenvalue problems as systems of
%   nonlinear equations.
%
%   The problem is, given a symmetric matrix A, to find an orthogonal
%   matrix Q and diagonal matrix D such that A = Q(T) D Q. %
%   The dimension n must be such that 0.5*(sqrt(1+4*n)-1) is an integer.
%
%   Example A: a diagonal matrix.
%
%   Source:  An idea by Nick Gould
%
%   Ph. Toint 21 VII 2021.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'eigenals';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 2 || round(sqrt(1+4*n)) ~= sqrt(1+4*n) )
         disp( [ ' ERROR in eigenals: n = ', int2str(n), ' is < 2!' ] )
      end
   else
      n = 110;
      p = 10;
   end
   p   = 0.5*(sqrt(1+4*n)-1);
   if ( abs( round( p ) - p ) > 1.e-15 )
      disp( [ ' ERROR in eigenals: n = ', int2str(n), ' is not such that 0.5*(sqrt(1+4*n)-1) is integer!' ] )
   end
   Q = speye(p);
   varargout{1} = full([ Q(:); ones(p,1) ]);      % x0
   varargout{2} = 0;                              % fstar
   varargout{3} = '';                             % xtype
   varargout{4} = [];                             % xlower
   varargout{5} = [];                             % xupper
   varargout{6} = [];                             % clower
   varargout{7} = [];                             % cupper
   varargout{8} = 'SUR2-AN-V-0';                  % class

case 'cpsstr'

   n   = varargin{1};
   p   = round(0.5*(sqrt(1+4*n)-1));
   ntr = p*(p+1)/2;
   eldom = cell( 1, 2*ntr );
   itr = 0;
   for j = 1:p
      for i = 1:j
        itr = itr + 1;
	if ( i == j )
           eldom{itr}     = [ [(i-1)*p+1:i*p] [p^2+1:n] ];
  	   eldom{ntr+itr} = [ (i-1)*p+1:i*p ];
	else
           eldom{itr}     = [ [(i-1)*p+1:i*p] [p^2+1:n] [(j-1)*p+1:j*p] ];
  	   eldom{ntr+itr} = [ [(i-1)*p+1:i*p] [(j-1)*p+1:j*p] ];
	end
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

   iel  = varargin{1};
   x    = varargin{2};
   n    = varargin{3};
   m    = length(x);
   p    = round(0.5*(sqrt(1+4*n)-1));
   twop = p+p;
   Ac   = zeros( p*(p+1)/2, 1 );
   Ic   = zeros( p*(p+1)/2, 1 );
   ij   = 1;
   for k = 1:p
      Ac( ij ) = k;
      Ic( ij ) = 1;
      ij       = ij + k + 1;
   end
   ntr = p*(p+1)/2;
   if ( iel <= ntr )
      if ( m == twop )
         riel = -Ac(iel);
         for k = 1:p
            riel = riel + x(k)^2*x(p+k);
         end
         varargout{1} = riel^2;
         if ( nargout > 1 )
            Jiel = zeros( m, 1 );
            for k = 1:p
	       Jiel(k)   = Jiel(k)   + 2*x(k)*x(p+k);
	       Jiel(p+k) = Jiel(p+k) + x(k)^2;
	    end
            varargout{2} = 2 * Jiel * riel;
	    if ( nargout > 2 )
	       Hiel = zeros( m, m );
               for k = 1:p
	          Hiel(k,k)   = Hiel(k,k)   + 2*x(p+k);
	          Hiel(k,p+k) = Hiel(k,p+k) + 2*x(k);
	          Hiel(p+k,k) = Hiel(k,p+k);
	       end
               varargout{3} = 2 * ( Jiel*Jiel.' + riel*Hiel );
	    end
	 end
      else
         riel = -Ac(iel);
         for k = 1:p
            riel = riel + x(k)*x(p+k)*x(twop+k);
         end
         varargout{1} = riel^2;
         if ( nargout > 1 )
            Jiel = zeros( m, 1 );
            for k = 1:p
	       Jiel(k)      = Jiel(k)      + x(p+k)*x(twop+k);
	       Jiel(p+k)    = Jiel(p+k)    +   x(k)*x(twop+k);
	       Jiel(twop+k) = Jiel(twop+k) +    x(k)*x(p+k);
	    end
            varargout{2} = 2 * Jiel * riel;
	    if ( nargout > 2 )
	       Hiel = zeros( m, m );
               for k = 1:p
	          Hiel(k,p+k)      = Hiel(k,p+k)      + x(twop+k);
	          Hiel(k,twop+k)   = Hiel(k,twop+k)   + x(p+k);
	          Hiel(p+k,k)      = Hiel(k,p+k);
	          Hiel(p+k,twop+k) = Hiel(p+k,twop+k) + x(k);
	          Hiel(twop+k,k)   = Hiel(k,twop+k);
	          Hiel(twop+k,p+k) = Hiel(p+k,twop+k);
	       end
               varargout{3} = 2 * ( Jiel*Jiel.' + riel*Hiel );
	    end
	 end
      end
   else
      riel = -Ic(iel-ntr);
      if ( m == p )
         for k = 1:p
            riel = riel + x(k)^2;
         end
         varargout{1} = riel^2;
         if ( nargout > 1 )
            Jiel = zeros( m, 1 );
            for k = 1:p
	        Jiel(k)   = Jiel(k) + 2*x(k);
	    end
            varargout{2} = 2 * Jiel * riel;
	    if ( nargout > 2 )
	       Hiel = zeros( m, m );
               for k = 1:p
	          Hiel(k,k) = Hiel(k,k) + 2;
	       end
               varargout{3} = 2 * ( Jiel*Jiel.' + riel*Hiel );
            end
	 end
      else
         for k = 1:p
            riel = riel + x(k)*x(p+k);
         end
         varargout{1} = riel^2;
         if ( nargout > 1 )
            Jiel = zeros( m, 1 );
            for k = 1:p
	        Jiel(k)   = Jiel(k)   + x(p+k);
	        Jiel(p+k) = Jiel(p+k) + x(k);
	    end
            varargout{2} = 2 * Jiel * riel;
	    if ( nargout > 2 )
	       Hiel = zeros( m, m );
               for k = 1:p
	          Hiel(k,p+k) = Hiel(k,p+k) + 1;
	          Hiel(p+k,k) = Hiel(k,p+k);
	       end
               varargout{3} = 2 * ( Jiel*Jiel.' + riel*Hiel );
            end
	 end
      end
      
   end
end

return

end