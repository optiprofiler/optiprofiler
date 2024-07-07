%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = crglvy( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The extended Cragg and Levy problem. This problem is a sum of m  sets
%   of 5 groups. The Hessian matrix is 7-diagonal.
%
%   The dimension n must satisfy n = 2*m+2 for some positive integer m.
%
%   Source:  problem 32 in
%      Ph. L. Toint,
%      "Test problems for partially separable optimization and results
%      for the routine PSPMIN",
%      Report 83/4, Department of Mathematics, FUNDP (Namur, B), 1983.
%   Also problem 18 in 
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   If the dimension is unspecified, the default n = 10 is chosen.
%
%   Ph. Toint 22 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'crglvy';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 4 || round( (n-2)/2 ) ~= (n-2)/2 )
         disp( [ ' ERROR in crglvy: n = ', int2str(n),' does not satisfy n = 2*m+2!' ] )
      end
   else
      n = 10;
   end
   varargout{1} = [ 1; 2*ones( n-1, 1) ];       % x0
   varargout{2} = [ 0 1.886566 1.5372D+01 ];    % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AY-V-0';                % class

case 'cpsstr'

   n   = varargin{1};
   m   = ( n - 2 ) / 2;
   iel = 1;
   eldom = cell( m+4, 1 );
   for i = 1:m
      eldom{ iel }   = [ 2*i-1   2*i ];
      eldom{ iel+1 } = [ 2*i   2*i+1 ];
      eldom{ iel+2 } = [ 2*i+1 2*i+2 ];
      eldom{ iel+3 } = [ 2*i-1 ];
      eldom{ iel+4 } = [ 2*i+2 ];
      iel = iel + 5;
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
   switch ( mod( iel, 5 ) )
   case 1
      e1 = exp(x(1));
      t  = e1 - x(2);
      Jt = [ e1; -1 ];
      Ht = [ e1, 0;
             0 , 0 ];
      varargout{1} = t^4;
      if ( nargout > 1 )
         varargout{2} = 4 * t^3 * Jt;
	 if ( nargout > 2 )
	    varargout{3} = 4*(3*t^2*Jt*Jt.'+t^3*Ht);
	 end
      end
   case 2
      t =  x(1) - x(2);
      Jt = [ 1; -1 ];
      varargout{1} = 100*t^6;
      if ( nargout > 1 )
         varargout{2} = 600 * t^5 * Jt;
         if ( nargout > 2 )
	    varargout{3} = 3000*t^4*Jt*Jt.';
	 end
      end
   case 3
      t =  x(1) - x(2);
      Jt =  [ 1; -1 ];
      varargout{1} = tan(t)^4;
      if ( nargout > 1 )
         varargout{2} = 4 * tan(t)^3 * (1/cos(t))^2 * Jt;
         if ( nargout > 2 )
	    varargout{3} = 4 * ( 3*tan(t)^2*(1/cos(t))^4 + 2*(tan(t)/cos(t))^3*sin(t) ) * Jt*Jt.';
	 end
      end
   case 4
      varargout{1} = x(1)^8;
      if ( nargout > 1 )
         varargout{2} = 8*x(1)^7;
         if ( nargout > 2 )
	    varargout{3} = 56*x(1)^6;
	 end
      end
   case 0
      varargout{1} = (x(1) - 1)^2;
      if ( nargout > 1 )
         varargout{2} = 2*(x(1)-1);
         if ( nargout > 2 )
	    varargout{3} = 2;
	 end
      end
   end

end

return

end