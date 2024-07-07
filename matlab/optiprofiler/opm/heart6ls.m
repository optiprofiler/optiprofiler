%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = heart6ls( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Dipole model of the heart (6 x 6 version).
%   This is a nonlinear least-squares.
%
%   Source:
%      J. E. Dennis, Jr., D. M. Gay, P. A. Vu,
%      "A New Nonlinear Equations Test Problem".
%      Tech. Rep. 83-16, Dept. of Math. Sci., Rice Univ., Houston, TX
%      June 1983, revised May 1985.
%
%   Ph. Toint 27 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'heart6ls';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]
%                   a  c  t  u  v  w
   varargout{1} = [ 0, 0, 1, 1, 1, 1 ]';        % x0
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AY-6-0';                % class

case 'cpsstr'

   eldom{1} = [ 1:6 ];
   eldom{2} = [ 1:6 ];
   eldom{3} = [ 1:6 ];
   eldom{4} = [ 1:6 ];
   eldom{5} = [ 1:6 ];
   eldom{6} = [ 1:6 ];
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

   iel    = varargin{1};
   x      = varargin{2};
   sum_Mx =  -0.816;
   sum_My =  -0.017;
   sum_A  =  -1.826;
   sum_B  =  -0.754;
   sum_C  =  -4.839;
   sum_D  =  -3.259;
   sum_E  = -14.023;
   sum_F  =  15.467;
   switch ( iel )
   case 1
      iv1 = [ 3 1 ];
      iv2 = [ 4 1 ];
      iv3 = [ 5 2 ];
      iv4 = [ 6 2 ];
      switch ( nargout )
      case 1
         fe1 = heart6ls( '2prod', x(iv1) );
         fe2 = heart6ls( 'vpv'  , x(iv2), sum_Mx );
         fe3 = heart6ls( '2prod', x(iv3) );
         fe4 = heart6ls( 'vpv'  , x(iv4), sum_My );
	 riel = fe1+fe2-fe3-fe4-sum_A;
	 varargout{1} = riel^2;
      case 2
         [ fe1, ge1 ] = heart6ls( '2prod', x(iv1) );
         [ fe2, ge2 ] = heart6ls( 'vpv'  , x(iv2), sum_Mx  );
         [ fe3, ge3 ] = heart6ls( '2prod', x(iv3) );
         [ fe4, ge4 ] = heart6ls( 'vpv'  , x(iv4), sum_My );
	 riel = fe1+fe2-fe3-fe4-sum_A;
	 varargout{1} = riel^2;
	 Jiel         = zeros(6,1);
	 Jiel(iv1)    = Jiel(iv1) + ge1;
	 Jiel(iv2)    = Jiel(iv2) + ge2;
	 Jiel(iv3)    = Jiel(iv3) - ge3;
	 Jiel(iv4)    = Jiel(iv4) - ge4;
         varargout{2} = 2 * Jiel * riel;
      case 3
         [ fe1, ge1, He1 ] = heart6ls( '2prod', x(iv1) );
         [ fe2, ge2, He2 ] = heart6ls( 'vpv'  , x(iv2), sum_Mx  );
         [ fe3, ge3, He3 ] = heart6ls( '2prod', x(iv3) );
         [ fe4, ge4, He4 ] = heart6ls( 'vpv'  , x(iv4), sum_My );
	 riel = fe1+fe2-fe3-fe4-sum_A;
	 varargout{1} = riel^2;
	 Jiel         = zeros(6,1);
	 Jiel(iv1)    = Jiel(iv1) + ge1;
	 Jiel(iv2)    = Jiel(iv2) + ge2;
	 Jiel(iv3)    = Jiel(iv3) - ge3;
	 Jiel(iv4)    = Jiel(iv4) - ge4;
         varargout{2} = 2 * Jiel * riel;
         Hiel = sparse( 6, 6 );
	 Hiel(iv1,iv1) = Hiel(iv1,iv1) + He1;
	 Hiel(iv2,iv2) = Hiel(iv2,iv2) + He2;
	 Hiel(iv3,iv3) = Hiel(iv3,iv3) - He3;
	 Hiel(iv4,iv4) = Hiel(iv4,iv4) - He4;
         varargout{3} = 2* (Jiel*Jiel.'+riel*Hiel);
      end
   case 2
      iv1 = [ 5 1 ];
      iv2 = [ 6 1 ];
      iv3 = [ 3 2 ];
      iv4 = [ 4 2 ];
      switch ( nargout )
      case 1
         fe1 = heart6ls( '2prod', x(iv1) );
         fe2 = heart6ls( 'vpv'  , x(iv2), sum_Mx );
         fe3 = heart6ls( '2prod', x(iv3) );
         fe4 = heart6ls( 'vpv'  , x(iv4), sum_My );
	 riel = fe1+fe2-fe3-fe4-sum_B;
	 varargout{1} = riel^2;
      case 2
         [ fe1, ge1 ] = heart6ls( '2prod', x(iv1) );
         [ fe2, ge2 ] = heart6ls( 'vpv'  , x(iv2), sum_Mx  );
         [ fe3, ge3 ] = heart6ls( '2prod', x(iv3) );
         [ fe4, ge4 ] = heart6ls( 'vpv'  , x(iv4), sum_My );
	 riel = fe1+fe2-fe3-fe4-sum_B;
 	 varargout{1} = riel^2;
	 Jiel         = zeros(6,1);
	 Jiel(iv1)    = Jiel(iv1) + ge1;
	 Jiel(iv2)    = Jiel(iv2) + ge2;
	 Jiel(iv3)    = Jiel(iv3) - ge3;
	 Jiel(iv4)    = Jiel(iv4) - ge4;
         varargout{2} = 2 * Jiel * riel;
      case 3
         [ fe1, ge1, He1 ] = heart6ls( '2prod', x(iv1) );
         [ fe2, ge2, He2 ] = heart6ls( 'vpv'  , x(iv2), sum_Mx  );
         [ fe3, ge3, He3 ] = heart6ls( '2prod', x(iv3) );
         [ fe4, ge4, He4 ] = heart6ls( 'vpv'  , x(iv4), sum_My );
	 riel = fe1+fe2-fe3-fe4-sum_B;
	 varargout{1} = riel^2;
	 Jiel         = zeros(6,1);
	 Jiel(iv1)    = Jiel(iv1) + ge1;
	 Jiel(iv2)    = Jiel(iv2) + ge2;
	 Jiel(iv3)    = Jiel(iv3) - ge3;
	 Jiel(iv4)    = Jiel(iv4) - ge4;
         varargout{2} = 2 * Jiel * riel;
         Hiel = sparse( 6, 6 );
	 Hiel(iv1,iv1) = Hiel(iv1,iv1) + He1;
	 Hiel(iv2,iv2) = Hiel(iv2,iv2) + He2;
	 Hiel(iv3,iv3) = Hiel(iv3,iv3) - He3;
	 Hiel(iv4,iv4) = Hiel(iv4,iv4) - He4;
         varargout{3} = 2* (Jiel*Jiel.'+riel*Hiel);
      end
   case 3
      iv1 = [ 1 3 5 ];
      iv2 = [ 2 3 5 ];
      iv3 = [ 1 4 6 ];
      iv4 = [ 2 4 6 ];
      switch ( nargout )
      case 1
         fe1 = heart6ls( 'adfsq', x(iv1) );
         fe2 = heart6ls( '3prod', x(iv2) );
         fe3 = heart6ls( 'pdfsq', x(iv3), sum_Mx );
         fe4 = heart6ls( 'p3prd', x(iv4), sum_My );
	 riel = fe1-2*fe2+fe3-2*fe4-sum_C;
	 varargout{1} = riel^2;
      case 2
         [ fe1, ge1 ] = heart6ls( 'adfsq', x(iv1) );
         [ fe2, ge2 ] = heart6ls( '3prod', x(iv2) );
         [ fe3, ge3 ] = heart6ls( 'pdfsq', x(iv3), sum_Mx );
         [ fe4, ge4 ] = heart6ls( 'p3prd', x(iv4), sum_My );
	 riel = fe1+fe2-fe3-fe4-sum_C;
	 varargout{1} = riel^2;
	 Jiel         = zeros(6,1);
	 Jiel(iv1)    = Jiel(iv1) + ge1;
	 Jiel(iv2)    = Jiel(iv2) - 2*ge2;
	 Jiel(iv3)    = Jiel(iv3) + ge3;
	 Jiel(iv4)    = Jiel(iv4) - 2*ge4;
         varargout{2} = 2 * Jiel * riel;
      case 3
         [ fe1, ge1, He1 ] = heart6ls( 'adfsq', x(iv1) );
         [ fe2, ge2, He2 ] = heart6ls( '3prod', x(iv2) );
         [ fe3, ge3, He3 ] = heart6ls( 'pdfsq', x(iv3), sum_Mx );
         [ fe4, ge4, He4 ] = heart6ls( 'p3prd', x(iv4), sum_My );
	 riel = fe1-2*fe2+fe3-2*fe4-sum_C;
	 varargout{1} = riel^2;
	 Jiel         = zeros(6,1);
	 Jiel(iv1)    = Jiel(iv1) + ge1;
	 Jiel(iv2)    = Jiel(iv2) - 2*ge2;
	 Jiel(iv3)    = Jiel(iv3) + ge3;
	 Jiel(iv4)    = Jiel(iv4) - 2*ge4;
         varargout{2} = 2 * Jiel * riel;
         Hiel = sparse( 6, 6 );
	 Hiel(iv1,iv1) = Hiel(iv1,iv1) + He1;
	 Hiel(iv2,iv2) = Hiel(iv2,iv2) - 2*He2;
	 Hiel(iv3,iv3) = Hiel(iv3,iv3) + He3;
	 Hiel(iv4,iv4) = Hiel(iv4,iv4) - 2*He4;
         varargout{3} = 2* (Jiel*Jiel.'+riel*Hiel);
      end
   case 4
      iv1 = [ 2 3 5 ];
      iv2 = [ 1 3 5 ];
      iv3 = [ 2 4 6 ];
      iv4 = [ 1 4 5 ];
      switch ( nargout )
      case 1
         fe1 = heart6ls( 'adfsq', x(iv1) );
         fe2 = heart6ls( '3prod', x(iv2) );
         fe3 = heart6ls( 'pdfsq', x(iv3), sum_My );
         fe4 = heart6ls( 'p3prd', x(iv4), sum_Mx );
	 riel = fe1+2*fe2+fe3+2*fe4-sum_D;
	 varargout{1} = riel^2;
      case 2
         [ fe1, ge1 ] = heart6ls( 'adfsq', x(iv1) );
         [ fe2, ge2 ] = heart6ls( '3prod', x(iv2) );
         [ fe3, ge3 ] = heart6ls( 'pdfsq', x(iv3), sum_Mx );
         [ fe4, ge4 ] = heart6ls( 'p3prd', x(iv4), sum_My );
	 riel = fe1+2*fe2+fe3+2*fe4-sum_D;
	 varargout{1} = riel^2;
	 Jiel         = zeros(6,1);
	 Jiel(iv1)    = Jiel(iv1) + ge1;
	 Jiel(iv2)    = Jiel(iv2) + 2*ge2;
	 Jiel(iv3)    = Jiel(iv3) + ge3;
	 Jiel(iv4)    = Jiel(iv4) + 2*ge4;
         varargout{2} = 2 * Jiel * riel;
      case 3
         [ fe1, ge1, He1 ] = heart6ls( 'adfsq', x(iv1) );
         [ fe2, ge2, He2 ] = heart6ls( '3prod', x(iv2) );
         [ fe3, ge3, He3 ] = heart6ls( 'pdfsq', x(iv3), sum_Mx );
         [ fe4, ge4, He4 ] = heart6ls( 'p3prd', x(iv4), sum_My );
	 riel = fe1+2*fe2+fe3+2*fe4-sum_D;
	 varargout{1} = riel^2;
	 Jiel         = zeros(6,1);
	 Jiel(iv1)    = Jiel(iv1) + ge1;
	 Jiel(iv2)    = Jiel(iv2) + 2*ge2;
	 Jiel(iv3)    = Jiel(iv3) + ge3;
	 Jiel(iv4)    = Jiel(iv4) + 2*ge4;
         varargout{2} = 2 * Jiel * riel;
         Hiel = sparse( 6, 6 );
	 Hiel(iv1,iv1) = Hiel(iv1,iv1) + He1;
	 Hiel(iv2,iv2) = Hiel(iv2,iv2) + 2*He2;
	 Hiel(iv3,iv3) = Hiel(iv3,iv3) + He3;
	 Hiel(iv4,iv4) = Hiel(iv4,iv4) + 2*He4;
         varargout{3} = 2* (Jiel*Jiel.'+riel*Hiel);
      end
   case 5
      iv1 = [ 1 3 5 ];
      iv2 = [ 2 5 3 ];
      iv3 = [ 1 4 6 ];
      iv4 = [ 2 6 4 ];
      switch ( nargout )
      case 1
         fe1 = heart6ls( '3dprd', x(iv1) );
         fe2 = heart6ls( '3dprd', x(iv2) );
         fe3 = heart6ls( 'd3prd', x(iv3), sum_Mx );
         fe4 = heart6ls( 'd3prd', x(iv4), sum_My );
	 riel = fe1+fe2+fe3+fe4-sum_E;
	 varargout{1} = riel^2;
      case 2
         [ fe1, ge1 ] = heart6ls( '3dprd', x(iv1) );
         [ fe2, ge2 ] = heart6ls( '3dprd', x(iv2) );
         [ fe3, ge3 ] = heart6ls( 'd3prd', x(iv3), sum_Mx );
         [ fe4, ge4 ] = heart6ls( 'd3prd', x(iv4), sum_My );
	 riel = fe1+fe2+fe3+fe4-sum_E;
	 varargout{1} = riel^2;
	 Jiel         = zeros(6,1);
	 Jiel(iv1)    = Jiel(iv1) + ge1;
	 Jiel(iv2)    = Jiel(iv2) + ge2;
	 Jiel(iv3)    = Jiel(iv3) + ge3;
	 Jiel(iv4)    = Jiel(iv4) + ge4;
         varargout{2} = 2 * Jiel * riel;
      case 3
         [ fe1, ge1, He1 ] = heart6ls( '3dprd', x(iv1) );
         [ fe2, ge2, He2 ] = heart6ls( '3dprd', x(iv2) );
         [ fe3, ge3, He3 ] = heart6ls( 'd3prd', x(iv3), sum_Mx );
         [ fe4, ge4, He4 ] = heart6ls( 'd3prd', x(iv4), sum_My );
	 riel = fe1+fe2+fe3+fe4-sum_E;
	 varargout{1} = riel^2;
	 Jiel         = zeros(6,1);
	 Jiel(iv1)    = Jiel(iv1) + ge1;
	 Jiel(iv2)    = Jiel(iv2) + ge2;
	 Jiel(iv3)    = Jiel(iv3) + ge3;
	 Jiel(iv4)    = Jiel(iv4) + ge4;
         varargout{2} = 2 * Jiel * riel;
         Hiel = sparse( 6, 6 );
	 Hiel(iv1,iv1) = Hiel(iv1,iv1) + He1;
	 Hiel(iv2,iv2) = Hiel(iv2,iv2) + He2;
	 Hiel(iv3,iv3) = Hiel(iv3,iv3) + He3;
	 Hiel(iv4,iv4) = Hiel(iv4,iv4) + He4;
         varargout{3} = 2* (Jiel*Jiel.'+riel*Hiel);
      end
   case 6
      iv1 = [ 2 3 5 ];
      iv2 = [ 1 3 5 ];
      iv3 = [ 2 4 6 ];
      iv4 = [ 1 6 4 ];
      switch ( nargout )
      case 1
         fe1 = heart6ls( '3dprd', x(iv1) );
         fe2 = heart6ls( '3dprd', x(iv2) );
         fe3 = heart6ls( 'd3prd', x(iv3), sum_My );
         fe4 = heart6ls( 'd3prd', x(iv4), sum_Mx );
	 riel = fe1-fe2+fe3-fe4-sum_F;
	 varargout{1} = riel^2;
      case 2
         [ fe1, ge1 ] = heart6ls( '3dprd', x(iv1) );
         [ fe2, ge2 ] = heart6ls( '3dprd', x(iv2) );
         [ fe3, ge3 ] = heart6ls( 'd3prd', x(iv3), sum_My );
         [ fe4, ge4 ] = heart6ls( 'd3prd', x(iv4), sum_Mx );
	 riel = fe1-fe2+fe3-fe4-sum_F;
	 varargout{1} = riel^2;
	 Jiel         = zeros(6,1);
	 Jiel(iv1)    = Jiel(iv1) + ge1;
	 Jiel(iv2)    = Jiel(iv2) - ge2;
	 Jiel(iv3)    = Jiel(iv3) + ge3;
	 Jiel(iv4)    = Jiel(iv4) - ge4;
         varargout{2} = 2 * Jiel * riel;
      case 3
         [ fe1, ge1, He1 ] = heart6ls( '3dprd', x(iv1) );
         [ fe2, ge2, He2 ] = heart6ls( '3dprd', x(iv2) );
         [ fe3, ge3, He3 ] = heart6ls( 'd3prd', x(iv3), sum_My );
         [ fe4, ge4, He4 ] = heart6ls( 'd3prd', x(iv4), sum_Mx );
	 riel = fe1-fe2+fe3-fe4-sum_F;
	 varargout{1} = riel^2;
	 Jiel         = zeros(6,1);
	 Jiel(iv1)    = Jiel(iv1) + ge1;
	 Jiel(iv2)    = Jiel(iv2) - ge2;
	 Jiel(iv3)    = Jiel(iv3) + ge3;
	 Jiel(iv4)    = Jiel(iv4) - ge4;
         varargout{2} = 2 * Jiel * riel;
         Hiel = sparse( 6, 6 );
	 Hiel(iv1,iv1) = Hiel(iv1,iv1) + He1;
	 Hiel(iv2,iv2) = Hiel(iv2,iv2) - He2;
	 Hiel(iv3,iv3) = Hiel(iv3,iv3) + He3;
	 Hiel(iv4,iv4) = Hiel(iv4,iv4) - He4;
         varargout{3} = 2* (Jiel*Jiel.'+riel*Hiel);
      end
   end

case '2prod' % varargout = [ fiel, giel, Hiel ]

   x = varargin{1};
   varargout{1} = x(1)*x(2);
   if ( nargout > 1 )
      varargout{2} = [ x(2); x(1) ];
      if ( nargout > 2 )
         varargout{3} = [ 0, 1 ;
	                  1, 0 ];
      end
   end
   
case '3prod' % varargout = [ fiel, giel, Hiel ]

   x = varargin{1};
   varargout{1} = x(1)*x(2)*x(3);
   if ( nargout > 1 )
      varargout{2} = [ x(2)*x(3); x(1)*x(3); x(1)*x(2) ];
      if ( nargout > 2 )
         varargout{3} = [  0  , x(3), x(2) ;
	                  x(3),   0 , x(1) ;
			  x(2), x(1),   0  ];
      end
   end
   
case 'vpv' % varargout = [ fiel, giel, Hiel ]

   x     = varargin{1};
   alpha = varargin{2};
   diff  = alpha - x(2);
   varargout{1} = x(1)*diff;
   if ( nargout > 1 )
      varargout{2} = [ diff; -x(1) ];
      if ( nargout > 2 )
         varargout{3} = [  0  , -1 ;
	                  -1  ,  0 ];
      end
   end

case 'adfsq' % varargout = [ fiel, giel, Hiel ]

   x     = varargin{1};
   dfsq  = x(2)^2-x(3)^2;
   varargout{1} = x(1)*dfsq;
   if ( nargout > 1 )
      twox = 2*x(1);
      varargout{2} = [ dfsq; twox*x(2); -twox*x(3)];
      if ( nargout > 2 )
         varargout{3} = [   0   ,  2*x(2), -2*x(3);
	                  2*x(2),   twox ,   0    ;
	                 -2*x(3),    0   , -twox   ];
      end
   end

case 'pdfsq' % varargout = [ fiel, giel, Hiel ]

   x     = varargin{1};
   alpha = varargin{2};
   diff  = alpha - x(1);
   dfsq  = x(2)^2-x(3)^2;
   twod  = 2*diff;
   varargout{1} = diff*dfsq;
   if ( nargout > 1 )
      varargout{2} = [ -dfsq; twod*x(2); -twod*x(3) ];
      if ( nargout > 2 )
         varargout{3} = [     0     , -2*x(2), 2*x(3);
	                  -2*x(2),     twod  ,     0    ;
			   2*x(3),       0   ,   -twod  ];
      end
   end

case 'p3prd'

   x     = varargin{1};
   alpha = varargin{2};
   diff  = alpha - x(1);
   varargout{1} = diff*x(2)*x(3);
   if ( nargout > 1 )
      varargout{2} = [ -x(2)*x(3); diff*x(3); diff*x(2) ];
      if ( nargout > 2 )
         varargout{3} = [     0, -x(3), -x(2);
	                  -x(3),   0  ,  diff;
			  -x(2), diff ,   0  ];
      end
   end
   
case '3dprd'

   x     = varargin{1};
   diff  = x(2)^2-3*x(3)^2;
   varargout{1} = diff*x(1)*x(2);
   if ( nargout > 1 )
      varargout{2} = [ x(2)*diff; diff*x(1)+2*x(1)*x(2)^2; -6*x(1)*x(2)*x(3) ];
      if ( nargout > 2 )
         varargout{3} = [     0       , diff+2*x(2)^2, -6*x(2)*x(3);
	                 diff+2*x(2)^2,  6*x(1)*x(2) , -6*x(1)*x(3);
			 -6*x(2)*x(3) , -6*x(1)*x(3) , -6*x(1)*x(2) ];
      end
   end
   
case 'd3prd'

   x     = varargin{1};
   alpha = varargin{2};
   dfsq  = x(2)^2-3*x(3)^2;
   diff  = alpha-x(1);
   varargout{1} = diff*x(2)*dfsq;
   if ( nargout > 1 )
      varargout{2} = [ -x(2)*dfsq; diff*(dfsq+2*x(2)^2); -6*x(2)*x(3)*diff ];
      if ( nargout > 2 )
         varargout{3} = [      0       , -dfsq-2*x(2)^2,  6*x(2)*x(3);
	                 -dfsq-2*x(2)^2,   6*x(2)*diff , -6*x(3)*diff;
			  6*x(2)*x(3)  ,  -6*x(3)*diff , -6*x(2)*diff ];
      end
   end
   
end

return

end