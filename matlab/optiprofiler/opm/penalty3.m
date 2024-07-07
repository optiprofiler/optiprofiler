%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = penalty3( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Third penalty problem.
%
%   Source: problem 71 (p. 81) in
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 26 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'penalty3';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 3 || round( n/2 ) ~= n/2 )
         disp( [ ' ERROR in penalty3: n = ', int2str(n), ' must be > 2 and even!' ] )
      end
   else
      n = 10;
   end
   varargout{1} = zeros( n , 1 );               % x0
   varargout{2} = 'unknown';                    % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'OUR2-AY-V-0';                % class

case 'cpsstr'

   n            = varargin{1};
   eldom{ 1 }   = [ 1:n ];
   eldom{ 2 }   = [ 1:n ];
   eldom{ 3 }   = [ 1:n ];
   eldom{ 4 }   = [ 1:n ];
   eldom{ 5 }   = [ 1:n/2 ];
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
   a   = 0.001;
   switch ( iel );
   case 1
      r  = [ x(1:n-2)+2*x(2:n-1)+10*x(3:n)-ones(n-2,1) ];
      R  = sum( r.^2 );
      en = exp(x(n));
      varargout{1} = a*(1+R*en);
      if ( nargout > 1 )
         Jr = zeros( n, n-2 );
         for i = 1:n-2
            Jr(1:n,i) = [zeros(i-1,1);1;2;10;zeros(n-i-2,1)];
	 end
         gR = 2*Jr*r;
         varargout{2} = a*(gR*en+R*[zeros(n-1,1);en]);
         if ( nargout > 2 )
	    varargout{3} = 2*a*en*Jr*Jr.';
	    varargout{3}(n,1:n) = varargout{3}(n,1:n)+a*gR.'*en;
	    varargout{3}(1:n,n) = varargout{3}(1:n,n)+a*gR*en;
	    varargout{3}(n,n)   = varargout{3}(n,n)  +a*R*en;
	 end
      end
   case 2
      s    = [ 2*x(1:n-2)+x(2:n-1)-3*ones(n-2,1) ];
      S    = sum( s.^2 );
      enm1 = exp(x(n-1));
      varargout{1} = a*S*enm1;
      if ( nargout > 1 )
         Js = zeros( n, n-2 );
         for i = 1:n-2
            Js(1:n,i)= [zeros(i-1,1);2;1;zeros(n-i-1,1)];
	 end
         gS = 2*Js*s;
         varargout{2} = a*(gS*enm1+S*[zeros(n-2,1);enm1;0]);
         if ( nargout > 2 )
	    varargout{3} = 2*a*enm1*Js*Js.';
	    varargout{3}(n-1,1:n) = varargout{3}(n-1,1:n)+a*gS.'*enm1;
	    varargout{3}(1:n,n-1) = varargout{3}(1:n,n-1)+a*gS*enm1;
	    varargout{3}(n-1,n-1) = varargout{3}(n-1,n-1)+a*S*enm1;
	 end
      end
   case 3
      r  = [ x(1:n-2)+2*x(2:n-1)+10*x(3:n)-ones(n-2,1) ];
      R  = sum( r.^2 );
      s    = [ 2*x(1:n-2)+x(2:n-1)-3*ones(n-2,1) ];
      S    = sum( s.^2 );
      varargout{1} = a*R*S;
      if ( nargout > 1 )
         Jr = zeros( n, n-2 );
         for i = 1:n-2
            Jr(1:n,i)= [zeros(i-1,1);1;2;10;zeros(n-i-2,1)];
	 end
         gR = 2*Jr*r;
         Js = zeros( n, n-2 );
         for i = 1:n-2
            Js(1:n,i)= [zeros(i-1,1);2;1;zeros(n-i-1,1)];
	 end
         gS = 2*Js*s;
         varargout{2} = a*(gR*S+R*gS);
	 if ( nargout > 2 )
	    HR = 2*Jr*Jr.';
	    HS = 2*Js*Js.';
	    varargout{3} = a * ( gR*gS.'+ HR*S + gS*gR.' + HS*R );
         end
      end
   case 4
      v = x.^2-n*ones(n,1);
      varargout{1} = sum( v.^2 );
      if ( nargout > 1 )
         Jv = 2*diag(x);
	 varargout{2} = 2*Jv*v;
         if ( nargout > 2 )
	    varargout{3} = 2 * ( Jv*Jv.'+ 2*diag(v) );
	 end
      end

   case 5
      w = x - ones(n/2,1);
      varargout{1} = sum( w.^2 );
      if ( nargout > 1 )
	 varargout{2} = 2 * w;
         if ( nargout > 2 )
	    varargout{3} = 2 * eye(n/2);
	 end
      end
   end

end

return

end