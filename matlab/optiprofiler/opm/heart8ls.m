%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = heart8ls( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Dipole model of the heart (8 x 8 version).
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

pname   = 'heart8ls';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]
%                   a  b  c  d  t  u  v  w
   varargout{1} = [ 0, 1, 0, 1, 1, 1, 1, 1 ]';  % x0
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-MN-8-0';                % class

case 'cpsstr'

   eldom{1} = [ 1 2 ];
   eldom{2} = [ 3 4 ];
   eldom{3} = [ 1:8 ];
   eldom{4} = [ 1:8 ];
   eldom{5} = [ 1:8 ];
   eldom{6} = [ 1:8 ];
   eldom{7} = [ 1:8 ];
   eldom{8} = [ 1:8 ];
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
   sum_Mx =  -0.69;
   sum_My =  -0.044;
   sum_A  =  -1.57;
   sum_B  =  -1.31;
   sum_C  =  -2.65;
   sum_D  =   2.0;
   sum_E  = -12.6;
   sum_F  =   9.48;
   switch ( iel )
   case 1
      riel = x(1)+x(2)-sum_Mx;
      varargout{1} = riel^2;
      if ( nargout > 1 )
	 Jiel  =  [ 1; 1 ];
         varargout{2} = 2 * Jiel * riel;
         if ( nargout > 2 )
            varargout{3} = 2 * Jiel*Jiel';
         end
      end
   case 2
      riel = x(1)+x(2)-sum_My;
      varargout{1} = riel^2;
      if ( nargout > 1 )
	 Jiel = [ 1; 1 ];
         varargout{2} = 2 * Jiel * riel;
         if ( nargout > 2 )
            varargout{3} = 2 * Jiel*Jiel';
         end
      end
   case 3
      iv1 = [ 5 1 ];
      iv2 = [ 6 2 ];
      iv3 = [ 7 3 ];
      iv4 = [ 8 4 ];
      switch ( nargout )
      case 1
         fe1  = heart8ls( '2prod', x(iv1) );
         fe2  = heart8ls( '2prod', x(iv2) );
         fe3  = heart8ls( '2prod', x(iv3) );
         fe4  = heart8ls( '2prod', x(iv4) );
	 riel = fe1+fe2-fe3-fe4 - sum_A;
	 varargout{1} = riel^2;
      case 2
         [ fe1, ge1 ] = heart8ls( '2prod', x(iv1) );
         [ fe2, ge2 ] = heart8ls( '2prod', x(iv2) );
         [ fe3, ge3 ] = heart8ls( '2prod', x(iv3) );
         [ fe4, ge4 ] = heart8ls( '2prod', x(iv4) );
	 riel         = fe1+fe2-fe3-fe4 - sum_A;
	 varargout{1} = riel^2;
	 Jiel         = zeros(8,1);
	 Jiel(iv1)    = Jiel(iv1) + ge1;
	 Jiel(iv2)    = Jiel(iv2) + ge2;
	 Jiel(iv3)    = Jiel(iv3) - ge3;
	 Jiel(iv4)    = Jiel(iv4) - ge4;
         varargout{2} = 2 * Jiel * riel;
      case 3
         [ fe1, ge1, He1 ] = heart8ls( '2prod', x(iv1) );
         [ fe2, ge2, He2 ] = heart8ls( '2prod', x(iv2) );
         [ fe3, ge3, He3 ] = heart8ls( '2prod', x(iv3) );
         [ fe4, ge4, He4 ] = heart8ls( '2prod', x(iv4) );
	 riel              = fe1+fe2-fe3-fe4 - sum_A;
	 varargout{1} = riel^2;
	 Jiel         = zeros(8,1);
	 Jiel(iv1)    = Jiel(iv1) + ge1;
	 Jiel(iv2)    = Jiel(iv2) + ge2;
	 Jiel(iv3)    = Jiel(iv3) - ge3;
	 Jiel(iv4)    = Jiel(iv4) - ge4;
         varargout{2} = 2 * Jiel * riel;
         Hiel = sparse( 8, 8 );
	 Hiel(iv1,iv1) = Hiel(iv1,iv1) + He1;
	 Hiel(iv2,iv2) = Hiel(iv2,iv2) + He2;
	 Hiel(iv3,iv3) = Hiel(iv3,iv3) - He3;
	 Hiel(iv4,iv4) = Hiel(iv4,iv4) - He4;
         varargout{3} = 2* (Jiel*Jiel.'+riel*Hiel);
      end
   case 4
      iv1 = [ 1 7 ];
      iv2 = [ 2 8 ];
      iv3 = [ 3 5 ];
      iv4 = [ 4 6 ];
      switch ( nargout )
      case 1
         fe1 = heart8ls( '2prod', x(iv1) );
         fe2 = heart8ls( '2prod', x(iv2) );
         fe3 = heart8ls( '2prod', x(iv3) );
         fe4 = heart8ls( '2prod', x(iv4) );
	 riel = fe1+fe2+fe3+fe4-sum_B;
	 varargout{1} = riel^2;
      case 2
         [ fe1, ge1 ] = heart8ls( '2prod', x(iv1) );
         [ fe2, ge2 ] = heart8ls( '2prod', x(iv2) );
         [ fe3, ge3 ] = heart8ls( '2prod', x(iv3) );
         [ fe4, ge4 ] = heart8ls( '2prod', x(iv4) );
	 riel = fe1+fe2+fe3+fe4-sum_B;
	 varargout{1} = riel^2;
	 Jiel         = zeros(8,1);
	 Jiel(iv1)    = Jiel(iv1) + ge1;
	 Jiel(iv2)    = Jiel(iv2) + ge2;
	 Jiel(iv3)    = Jiel(iv3) + ge3;
	 Jiel(iv4)    = Jiel(iv4) + ge4;
         varargout{2} = 2 * Jiel * riel;
      case 3
         [ fe1, ge1, He1 ] = heart8ls( '2prod', x(iv1) );
         [ fe2, ge2, He2 ] = heart8ls( '2prod', x(iv2) );
         [ fe3, ge3, He3 ] = heart8ls( '2prod', x(iv3) );
         [ fe4, ge4, He4 ] = heart8ls( '2prod', x(iv4) );
	 riel = fe1+fe2+fe3+fe4-sum_B;
	 varargout{1} = riel^2;
	 Jiel         = zeros(8,1);
	 Jiel(iv1)    = Jiel(iv1) + ge1;
	 Jiel(iv2)    = Jiel(iv2) + ge2;
	 Jiel(iv3)    = Jiel(iv3) + ge3;
	 Jiel(iv4)    = Jiel(iv4) + ge4;
         varargout{2} = 2 * Jiel * riel;
         Hiel = sparse( 8, 8 );
	 Hiel(iv1,iv1) = Hiel(iv1,iv1) + He1;
	 Hiel(iv2,iv2) = Hiel(iv2,iv2) + He2;
	 Hiel(iv3,iv3) = Hiel(iv3,iv3) + He3;
	 Hiel(iv4,iv4) = Hiel(iv4,iv4) + He4;
         varargout{3} = 2* (Jiel*Jiel.'+riel*Hiel);
      end
   case 5
      iv1 = [ 1 5 7 ];
      iv2 = [ 3 5 7 ];
      iv3 = [ 2 6 8 ];
      iv4 = [ 4 6 8 ];
      switch ( nargout )
      case 1
         fe1 = heart8ls( 'adfsq', x(iv1) );
         fe2 = heart8ls( '3prod', x(iv2) );
         fe3 = heart8ls( 'adfsq', x(iv3) );
         fe4 = heart8ls( '3prod', x(iv4) );
	 riel = fe1-2*fe2+fe3-2*fe4-sum_C;
	 varargout{1} = riel^2;
      case 2
         [ fe1, ge1 ] = heart8ls( 'adfsq', x(iv1) );
         [ fe2, ge2 ] = heart8ls( '3prod', x(iv2) );
         [ fe3, ge3 ] = heart8ls( 'adfsq', x(iv3) );
         [ fe4, ge4 ] = heart8ls( '3prod', x(iv4) );
	 riel = fe1-2*fe2+fe3-2*fe4-sum_C;
	 varargout{1} = riel^2;
	 Jiel         = zeros(8,1);
	 Jiel(iv1)    = Jiel(iv1) + ge1;
	 Jiel(iv2)    = Jiel(iv2) - 2*ge2;
	 Jiel(iv3)    = Jiel(iv3) + ge3;
	 Jiel(iv4)    = Jiel(iv4) - 2*ge4;
         varargout{2} = 2 * Jiel * riel;
      case 3
         [ fe1, ge1, He1 ] = heart8ls( 'adfsq', x(iv1) );
         [ fe2, ge2, He2 ] = heart8ls( '3prod', x(iv2) );
         [ fe3, ge3, He3 ] = heart8ls( 'adfsq', x(iv3) );
         [ fe4, ge4, He4 ] = heart8ls( '3prod', x(iv4) );
	 riel = fe1-2*fe2+fe3-2*fe4-sum_C;
	 varargout{1} = riel^2;
	 Jiel         = zeros(8,1);
	 Jiel(iv1)    = Jiel(iv1) + ge1;
	 Jiel(iv2)    = Jiel(iv2) - 2*ge2;
	 Jiel(iv3)    = Jiel(iv3) + ge3;
	 Jiel(iv4)    = Jiel(iv4) - 2*ge4;
         varargout{2} = 2 * Jiel * riel;
         Hiel = sparse( 8, 8 );
	 Hiel(iv1,iv1) = Hiel(iv1,iv1) + He1;
	 Hiel(iv2,iv2) = Hiel(iv2,iv2) - 2*He2;
	 Hiel(iv3,iv3) = Hiel(iv3,iv3) + He3;
	 Hiel(iv4,iv4) = Hiel(iv4,iv4) - 2*He4;
         varargout{3} = 2* (Jiel*Jiel.'+riel*Hiel);
      end
   case 6
      iv1 = [ 3 5 7 ];
      iv2 = [ 1 5 7 ];
      iv3 = [ 4 6 8 ];
      iv4 = [ 2 6 8 ];
      switch ( nargout )
      case 1
         fe1 = heart8ls( 'adfsq', x(iv1) );
         fe2 = heart8ls( '3prod', x(iv2) );
         fe3 = heart8ls( 'adfsq', x(iv3) );
         fe4 = heart8ls( '3prod', x(iv4) );
	 riel = fe1+2*fe2+fe3+2*fe4-sum_D;
	 varargout{1} = riel^2;
      case 2
         [ fe1, ge1 ] = heart8ls( 'adfsq', x(iv1) );
         [ fe2, ge2 ] = heart8ls( '3prod', x(iv2) );
         [ fe3, ge3 ] = heart8ls( 'adfsq', x(iv3) );
         [ fe4, ge4 ] = heart8ls( '3prod', x(iv4) );
	 riel = fe1+2*fe2+fe3+2*fe4-sum_D;
	 varargout{1} = riel^2;
	 Jiel         = zeros(8,1);
	 Jiel(iv1)    = Jiel(iv1) + ge1;
	 Jiel(iv2)    = Jiel(iv2) + 2*ge2;
	 Jiel(iv3)    = Jiel(iv3) + ge3;
	 Jiel(iv4)    = Jiel(iv4) + 2*ge4;
         varargout{2} = 2 * Jiel * riel;
      case 3
         [ fe1, ge1, He1 ] = heart8ls( 'adfsq', x(iv1) );
         [ fe2, ge2, He2 ] = heart8ls( '3prod', x(iv2) );
         [ fe3, ge3, He3 ] = heart8ls( 'adfsq', x(iv3) );
         [ fe4, ge4, He4 ] = heart8ls( '3prod', x(iv4) );
	 riel = fe1+2*fe2+fe3+2*fe4-sum_D;
	 varargout{1} = riel^2;
	 Jiel         = zeros(8,1);
	 Jiel(iv1)    = Jiel(iv1) + ge1;
	 Jiel(iv2)    = Jiel(iv2) + 2*ge2;
	 Jiel(iv3)    = Jiel(iv3) + ge3;
	 Jiel(iv4)    = Jiel(iv4) + 2*ge4;
         varargout{2} = 2 * Jiel * riel;
         Hiel = sparse( 8, 8 );
	 Hiel(iv1,iv1) = Hiel(iv1,iv1) + He1;
	 Hiel(iv2,iv2) = Hiel(iv2,iv2) + 2*He2;
	 Hiel(iv3,iv3) = Hiel(iv3,iv3) + He3;
	 Hiel(iv4,iv4) = Hiel(iv4,iv4) + 2*He4;
         varargout{3} = 2* (Jiel*Jiel.'+riel*Hiel);
      end
   case 7
      iv1 = [ 1 5 7 ];
      iv2 = [ 3 5 7 ];
      iv3 = [ 2 6 8 ];
      iv4 = [ 4 6 8 ];
      switch ( nargout )
      case 1
         fe1 = heart8ls( '3prod', x(iv1) );
         fe2 = heart8ls( '3prod', x(iv2) );
         fe3 = heart8ls( '3prod', x(iv3) );
         fe4 = heart8ls( '3prod', x(iv4) );
	 riel = fe1+fe2+fe3+fe4-sum_E;
	 varargout{1} = riel^2;
      case 2
         [ fe1, ge1 ] = heart8ls( '3prod', x(iv1) );
         [ fe2, ge2 ] = heart8ls( '3prod', x(iv2) );
         [ fe3, ge3 ] = heart8ls( '3prod', x(iv3) );
         [ fe4, ge4 ] = heart8ls( '3prod', x(iv4) );
	 riel = fe1+fe2+fe3+fe4-sum_E;
	 varargout{1} = riel^2;
	 Jiel         = zeros(8,1);
	 Jiel(iv1)    = Jiel(iv1) + ge1;
	 Jiel(iv2)    = Jiel(iv2) + ge2;
	 Jiel(iv3)    = Jiel(iv3) + ge3;
	 Jiel(iv4)    = Jiel(iv4) + ge4;
         varargout{2} = 2 * Jiel * riel;
      case 3
         [ fe1, ge1, He1 ] = heart8ls( '3prod', x(iv1) );
         [ fe2, ge2, He2 ] = heart8ls( '3prod', x(iv2) );
         [ fe3, ge3, He3 ] = heart8ls( '3prod', x(iv3) );
         [ fe4, ge4, He4 ] = heart8ls( '3prod', x(iv4) );
	 riel = fe1+fe2+fe3+fe4-sum_E;
	 varargout{1} = riel^2;
	 Jiel         = zeros(8,1);
	 Jiel(iv1)    = Jiel(iv1) + ge1;
	 Jiel(iv2)    = Jiel(iv2) + ge2;
	 Jiel(iv3)    = Jiel(iv3) + ge3;
	 Jiel(iv4)    = Jiel(iv4) + ge4;
         varargout{2} = 2 * Jiel * riel;
         Hiel = sparse( 8, 8 );
	 Hiel(iv1,iv1) = Hiel(iv1,iv1) + He1;
	 Hiel(iv2,iv2) = Hiel(iv2,iv2) + He2;
	 Hiel(iv3,iv3) = Hiel(iv3,iv3) + He3;
	 Hiel(iv4,iv4) = Hiel(iv4,iv4) + He4;
         varargout{3} = 2* (Jiel*Jiel.'+riel*Hiel);
      end
   case 8
      iv1 = [ 3 5 7 ];
      iv2 = [ 1 5 7 ];
      iv3 = [ 4 6 8 ];
      iv4 = [ 2 6 8 ];
      switch ( nargout )
      case 1
         fe1 = heart8ls( '3prod', x(iv1) );
         fe2 = heart8ls( '3prod', x(iv2) );
         fe3 = heart8ls( '3prod', x(iv3) );
         fe4 = heart8ls( '3prod', x(iv4) );
	 riel = fe1-fe2+fe3-fe4-sum_F;
	 varargout{1} = riel^2;
      case 2
         [ fe1, ge1 ] = heart8ls( '3prod', x(iv1) );
         [ fe2, ge2 ] = heart8ls( '3prod', x(iv2) );
         [ fe3, ge3 ] = heart8ls( '3prod', x(iv3) );
         [ fe4, ge4 ] = heart8ls( '3prod', x(iv4) );
	 riel = fe1-fe2+fe3-fe4-sum_F;
	 varargout{1} = riel^2;
	 Jiel         = zeros(8,1);
	 Jiel(iv1)    = Jiel(iv1) + ge1;
	 Jiel(iv2)    = Jiel(iv2) - ge2;
	 Jiel(iv3)    = Jiel(iv3) + ge3;
	 Jiel(iv4)    = Jiel(iv4) - ge4;
         varargout{2} = 2 * Jiel * riel;
      case 3
         [ fe1, ge1, He1 ] = heart8ls( '3prod', x(iv1) );
         [ fe2, ge2, He2 ] = heart8ls( '3prod', x(iv2) );
         [ fe3, ge3, He3 ] = heart8ls( '3prod', x(iv3) );
         [ fe4, ge4, He4 ] = heart8ls( '3prod', x(iv4) );
	 riel = fe1-fe2+fe3-fe4-sum_F;
	 varargout{1} = riel^2;
	 Jiel         = zeros(8,1);
	 Jiel(iv1)    = Jiel(iv1) + ge1;
	 Jiel(iv2)    = Jiel(iv2) - ge2;
	 Jiel(iv3)    = Jiel(iv3) + ge3;
	 Jiel(iv4)    = Jiel(iv4) - ge4;
         varargout{2} = 2 * Jiel * riel;
         Hiel = sparse( 8, 8 );
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

end

return

end