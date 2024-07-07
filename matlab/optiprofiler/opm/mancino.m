%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = mancino( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Mancino's function with variable dimension.
%
%   Source:
%   E. Spedicato,
%      "Computational experience with quasi-Newton algorithms for
%      minimization problems of moderate size",
%      Report N-175, CISE, Milano, 1975.
%   Also problem51 (p. 72) in
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   The definitions
%     s_{i,j} = sin log v_{i,j}   and s_{i,j} = cos log v_{i,j}
%   have been used.  It seems that the additional exponent alpha
%   in Buckley is a typo. The values alpha = 5, beta = 14 and gamma = 3
%   are used. A simple starting point is also chosen.

%   Also problem 31 (p. 70) in 
%
%   Ph. Toint 21 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname = 'mancino';

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   
   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 2 )
         disp( [ ' ERROR in mancino: n = ', int2str(n), ' must be > 1!' ] )
      end
   else
      n = 10;
   end
   varargout{1} = (1/n)*ones( n , 1 );          % x0
   varargout{2} = 'unknown';                    % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AN-V-0';                % class

case 'cpsstr'

   cpsstr.name  = pname;
   cpsstr.eldom = { [1:varargin{1}] };
   varargout{1} = cpsstr;
   
case 'objf'   % varargout = [ f, g, H ]

   x = varargin{1};
   n = length( x );
   varargout = opm_eval_cpsf( pname, 'elobjf', x, { [ 1:n ] }, nargout, n );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   x = varargin{2};
   n = varargin{3};
   varargout{1} = 0;
   if ( nargout > 1 )
      varargout{2} = zeros( n, 1 );
      if ( nargout > 2 )
         varargout{3} = zeros( n, n );
      end
   end
   for i = 1:n
      a     = x.^2 + i*[1:n]'.^(-1);
      vi    = sqrt( a );
      si    = sin( log( vi ) );
      ci    = cos( log( vi ) );
      hi    = vi.*( si + ci );
      ri    = sum( hi(1:i-1) )+ sum( hi(i+1:n) ) + 14 * n * x(i) + ( i - n/2 )^3;
      varargout{1} = varargout{1} + ri^2;
      if ( nargout > 1 )
         d1vi  =  x./vi;
         wi    =  d1vi./vi;
	 d1si  =  ci.*wi;
	 d1ci  = -si.*wi;
	 d1hi  = d1vi.*(si+ci)+vi.*(d1si+d1ci);
         Ji    = [ d1hi(1:i-1); 14*n; d1hi(i+1:n) ];
         varargout{2} = varargout{2} + 2* Ji * ri;
         if ( nargout > 2 )
	    d2vi =  -( (x./vi).^2 - 1)./vi;
	    d2si = -si.*wi.^2 - ci.*wi.^2 + ci.*d2vi./vi;
	    d2ci = -ci.*wi.^2 + si.*wi.^2 - si.*d2vi./vi;
	    d2hi = d2vi.*(si+ci)+2*d1vi.*(d1si+d1ci)+vi.*(d2si+d2ci);
	    Hi   = diag([d2hi(1:i-1);0;d2hi(i+1:n)]);
	    varargout{3} = varargout{3}+2*(Ji*Ji.'+ri*Hi);
         end
      end
   end
end

return

end