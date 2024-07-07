%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = hydc20ls( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The hydrocarbon-20 problem by Fletcher.
%   This is a least-squares version of problem HYDCAR20.
%
%   Source: Problem 2b in
%   J.J. More',"A collection of nonlinear model problems"
%   Proceedings of the AMS-SIAM Summer Seminar on the Computational
%   Solution of Nonlinear Systems of Equations, Colorado, 1988.
%   Argonne National Laboratory MCS-P60-0289, 1989.
%
%   SIF input : N. Gould and Ph. Toint, Feb 1991.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'hydc20ls';
problem = str2func( pname );
nn      = 20;
mm      = 3;
kk      = 9;

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n ~= 99 )
         disp( [ ' ERROR in hydc20ls: n = ', int2str(n), ' is not equal to 99!' ] )
      end
   end

%   Note: the correspondance between the SIF file variables and the present
%   ones is given by
%       X(I,J) (I=0:N-1,J=1:M)    =>   x(M*I+J)       (I=0:N-1, j=1:MM)
%       T(I)   (I=0:N-1)          =>   x(M*N+I+1)     (I=0:N-1)
%       V(I)   (I=0:N-2)          =>   x((M+1)*N+I+1) (I=0:N-2)

%   Starting point

%                 J = 1     2     3     I
   varargout{1} = [  0.0;  0.3;  0.1; % 0
                     0.0;  0.3;  0.9; % 1
		     0.01; 0.3;  0.9; % 2
                     0.02; 0.4;  0.8; % 3
		     0.05; 0.4;  0.8; % 4
		     0.07; 0.45; 0.8; % 5
		     0.09; 0.5;  0.7; % 6
		     0.1;  0.5;  0.7; % 7
		     0.15; 0.5;  0.6; % 8
		     0.2;  0.5;  0.6; % 9
		     0.25; 0.6;  0.5; % 10
		     0.3;  0.6;  0.5; % 11
		     0.35; 0.6;  0.5; % 12
                     0.4;  0.6;  0.4; % 13
		     0.4;  0.7;  0.4; % 14
		     0.42; 0.7;  0.3; % 15
		     0.45; 0.75; 0.3; % 16
		     0.45; 0.75; 0.2; % 17
		     0.5;  0.8;  0.1; % 18
		     0.5;  0.8;  0.0; % 19
		       100*ones(20,1);     %  T
		       300*ones(19,1)  ];  %  V
		       
   varargout{2} = 0;                              % fstar
   varargout{3} = '';                             % xtype
   varargout{4} = [];                             % xlower
   varargout{5} = [];                             % xupper
   varargout{6} = [];                             % clower
   varargout{7} = [];                             % cupper
   varargout{8} = 'SUR2-AN-99-0';                 % class

case 'cpsstr'

   iel = 1;
   for j = 1:mm     % 2.1-(j):  X(0,J) X(1,J) T(0)      V(0)
      eldom{iel} =               [ j    mm+j mm*nn+1 (mm+1)*nn+1 ];   
      iel = iel + 1;
   end
   for j = 1:mm     % 2.3-(j): X(N-2,J)      X(N-1,J)    T(N-2)
      eldom{iel} =           [ mm*(nn-2)+j, mm*(nn-1)+j, 4*nn-1 ];
      iel = iel + 1;
   end
   for j = 1:mm
      for i = 1:nn-2    % 2.2-(i,j): X(I-1,J)    X(I,J)  X(I+1,J)   T(I-1)  T(I)     V(I-1)   V(I)
         eldom{iel} =               [ mm*(i-1)+j mm*i+j mm*(i+1)+j mm*nn+i mm*nn+i+1 4*nn+i 4*nn+i+1 ];
         iel = iel + 1;
      end
   end
   for i = 0:nn-1       % 2.7-(i):   X(I,1) X(I,2) X(I,3)   T(I)
      eldom{iel} =                [ [mm*i+1:mm*i+mm] mm*nn+i+1 ];      
      iel = iel + 1;
   end
                        % 2.8 : X(0,[1:M]) X(1,[1:M])    T(0)    T(1)      V(0)
   eldom{iel} =                 [ [1:mm]  [mm+1:mm+mm] mm*nn+1 mm*nn+2 (mm+1)*nn+1];    
   iel = iel + 1;
   for i = 1:nn-2       % 2.9-(i):     X(I-1,J)                 X(I,J)              X(I+1,J)   
      eldom{iel} =              [  [mm*(i-1)+1:mm*(i-1)+mm]  [mm*i+1:mm*i+mm]   [mm*(i+1)+1:mm*(i+1)+mm] ];
                        %             T(I-1)     T(I)      T(i+1)     V(I-1)    V(I)
      eldom{iel} = [ eldom{iel}      mm*nn+i  mm*nn+i+1   mm*nn+i+2  4*nn+i   4*nn+i+1 ];       
      iel = iel + 1;
   end
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

%  Antoine constants, liquid and vapour enthalpies, etc

   a   = [  9.647    9.953     9.466  ];
   b   = [ -2998.00 -3448.10 -3347.25 ];
   c   = [   230.66   235.88   215.31 ];
   al  = [     0        0         0   ];
   al1 = [    37.6     48.2     45.4  ];
   al2 = [     0        0         0   ];
   be  = [  8425.0    9395.0   10466.0];
   be1 = [   24.2      35.6     31.9  ];
   be2 = [     0        0         0   ];
   fl  = [   30.0     30.0      40.0  ];
%  fv  = [     0        0         0   ];
   tf  = 100;
   bb  = 40;
   d   = 60;
   q   = 2500000.0;
   smallHF = 0;
   bigHF   = 0;
   for j = 1:mm
      smallHF = smallHF + fl(j)*( al(j)+al1(j)*tf+al2(j)*tf^2);
%     bigHF   = bigHF   + fv(j)*(be(j)+be1(j)*tf+ be2(j)*tf^2);
   end
   
   iel   = varargin{1};
   x     = varargin{2};
   niel  = length( x );
   if ( iel <= mm )                 % 2.1-(j)
      ie11 = [ 2 4 ];
      ie12 = [ 4 1 3 ];
      pe11 = [ -1 bb ];
      pe12 = [ -1 1 a(iel) b(iel) c(iel) ];
      switch( nargout )
      case 1
         fe11j = hydc20ls(    '2prod', x(ie11), pe11 );
         fe12j = hydc20ls( 'exp3prod', x(ie12), pe12 );
         riel  = bb*x(1) + fe11j + fe12j;
         varargout{1} = 0.0001*riel^2;
      case 2
         Jiel  = zeros( niel, 1 );
         [ fe11j, ge11j ] = hydc20ls(    '2prod', x(ie11), pe11 );
         [ fe12j, ge12j ] = hydc20ls( 'exp3prod', x(ie12), pe12 );
         riel         = bb*x(1) + fe11j + fe12j;
         varargout{1} = 0.0001*riel^2;
	 Jiel(1)      = Jiel(1) + bb;
	 Jiel(ie11)   = Jiel(ie11) + ge11j;
	 Jiel(ie12)   = Jiel(ie12) + ge12j;
	 varargout{2} = 0.0002*Jiel*riel;
      case 3
         Jiel    = zeros( niel, 1 );
	 Hiel    = zeros( niel, niel );
         [ fe11j, ge11j, He11j ] = hydc20ls(    '2prod', x(ie11), pe11 );
         [ fe12j, ge12j, He12j ] = hydc20ls( 'exp3prod', x(ie12), pe12 );
         riel            = bb*x(1) + fe11j + fe12j;
         varargout{1}    = 0.0001*riel^2;
	 Jiel(1)         = Jiel(1) + bb;
	 Jiel(ie11)      = Jiel(ie11) + ge11j;
	 Jiel(ie12)      = Jiel(ie12) + ge12j;
	 varargout{2}    = 0.0002*Jiel*riel;
	 Hiel(ie11,ie11) = Hiel(ie11,ie11) + He11j;
	 Hiel(ie12,ie12) = Hiel(ie12,ie12) + He12j;
	 varargout{3}    = 0.0002*(Jiel*Jiel'+riel*Hiel);
      end
   elseif ( iel <= mm+mm )
      j = iel-mm;
      ie31 = [ 1 3 ];
      pe31 = [ -1 1 a(j) b(j) c(j) ];
      switch( nargout )
      case 1
         fe31j = hydc20ls( 'exp2prod', x(ie31), pe31 );
	 riel  = -x(2) + fe31j;
	 varargout{1} = riel^2;
      case 2
         Jiel  = zeros( niel, 1 );
         [ fe31j, ge31j ] = hydc20ls( 'exp2prod', x(ie31), pe31 );
	 riel         = -x(2) + fe31j;
	 varargout{1} = riel^2;
	 Jiel(2)      = Jiel(2) - 1;
	 Jiel(ie31)   = Jiel(ie31) + ge31j;
	 varargout{2} = 2*Jiel*riel;
      case 3
         Jiel  = zeros( niel, 1 );
         Hiel  = zeros( niel, niel );
         [ fe31j, ge31j, He31j ] = hydc20ls( 'exp2prod', x(ie31), pe31 );
	 riel            = -x(2) + fe31j;
	 varargout{1}    = riel^2;
	 Jiel(2)         = Jiel(2) - 1;
	 Jiel(ie31)      = Jiel(ie31) + ge31j;
	 varargout{2}    = 2*Jiel*riel;
	 Hiel(ie31,ie31) = Hiel(ie31,ie31) + He31j;
	 varargout{3}    = 2*(Jiel*Jiel' + riel*Hiel);
      end
   elseif ( iel <= mm*nn )
%      i = floor( ( iel-2*mm)/mm );
%      if ( i == kk )
%         riel = - fl(j);
%      elseif ( i == kk+1 )
%         riel = - fv(j);
%      else
         riel = 0;
%      end
      j = mod( iel-2*mm, mm );
      if ( j == 0 )
         j = mm;
      end
      ie21 = [ 3 7 ];
      ie22 = [ 6 1 4 ];
      ie23 = [ 2 6 ];
      ie24 = [ 7 2 5 ];
      pe21 = [ -1 0 ];
      pe22 = [ -1 1 a(j) b(j) c(j) ];
      pe23 = [  1 0 ];
      pe24 = [  1 1 a(j) b(j) c(j) ];
      switch( nargout )
      case 1
         fe21j = hydc20ls(    '2prod', x(ie21), pe21 );
         fe22j = hydc20ls( 'exp3prod', x(ie22), pe22 );
         fe23j = hydc20ls(    '2prod', x(ie23), pe23 );
         fe24j = hydc20ls( 'exp3prod', x(ie24), pe24 );
	 riel  = riel + fe21j + fe22j + fe23j + fe24j;
	 varargout{1} = 0.0001*riel^2;
      case 2
         Jiel = zeros( niel, 1 );
         [ fe21j, ge21j ] = hydc20ls(    '2prod', x(ie21), pe21 );
         [ fe22j, ge22j ] = hydc20ls( 'exp3prod', x(ie22), pe22 );
         [ fe23j, ge23j ] = hydc20ls(    '2prod', x(ie23), pe23 );
         [ fe24j, ge24j ] = hydc20ls( 'exp3prod', x(ie24), pe24 );
	 riel         = riel + fe21j + fe22j + fe23j + fe24j;
	 varargout{1} = 0.0001*riel^2;
         Jiel(ie21)   = Jiel(ie21) + ge21j;
         Jiel(ie22)   = Jiel(ie22) + ge22j;
         Jiel(ie23)   = Jiel(ie23) + ge23j;
         Jiel(ie24)   = Jiel(ie24) + ge24j;
	 varargout{2} = 0.0002*Jiel*riel;
      case 3
         Jiel = zeros( niel, 1 );
	 Hiel = zeros( niel, niel );
         [ fe21j, ge21j, He21j ] = hydc20ls(    '2prod', x(ie21), pe21 );
         [ fe22j, ge22j, He22j ] = hydc20ls( 'exp3prod', x(ie22), pe22 );
         [ fe23j, ge23j, He23j ] = hydc20ls(    '2prod', x(ie23), pe23 );
         [ fe24j, ge24j, He24j ] = hydc20ls( 'exp3prod', x(ie24), pe24 );
	 riel         = riel + fe21j + fe22j + fe23j + fe24j;
	 varargout{1} = 0.0001*riel^2;
         Jiel(ie21)   = Jiel(ie21) + ge21j;
         Jiel(ie22)   = Jiel(ie22) + ge22j;
         Jiel(ie23)   = Jiel(ie23) + ge23j;
         Jiel(ie24)   = Jiel(ie24) + ge24j;
	 varargout{2} = 0.0002*Jiel*riel;
         Hiel(ie21,ie21) = Hiel(ie21,ie21) + He21j;
         Hiel(ie22,ie22) = Hiel(ie22,ie22) + He22j;
         Hiel(ie23,ie23) = Hiel(ie23,ie23) + He23j;
         Hiel(ie24,ie24) = Hiel(ie24,ie24) + He24j;
	 varargout{3}    = 0.0002*(Jiel*Jiel' + riel*Hiel);
      end
    elseif ( iel <= (mm+1)*nn )
%      i = iel-mm*nn;
      switch( nargout )
      case 1
         riel = -1;
	 for j = 1:mm
	    ie71 = [ j 4 ];
            pe71 = [ 1 1 a(j) b(j) c(j) ];
            fe71j = hydc20ls( 'exp2prod', x(ie71), pe71 );
	    riel  = riel + fe71j;
	 end
	 varargout{1} = riel^2;
      case 2
	 riel = -1;
         Jiel  = zeros( niel, 1 );
	 for j = 1:mm
	    ie71 = [ j 4 ];
            pe71 = [ 1 1 a(j) b(j) c(j) ];
            [ fe71j, ge71j ] = hydc20ls( 'exp2prod', x(ie71), pe71 );
	    riel       = riel + fe71j;
   	    Jiel(ie71) = Jiel(ie71) + ge71j;
	 end
	 varargout{1} = riel^2;
	 varargout{2} = 2*Jiel*riel;
      case 3
         riel = -1;
         Jiel  = zeros( niel, 1 );
         Hiel  = zeros( niel, niel );
         for j = 1:mm
	    ie71 = [ j 4 ];
            pe71 = [ 1 1 a(j) b(j) c(j) ];
            [ fe71j, ge71j, He71j ] = hydc20ls( 'exp2prod', x(ie71), pe71 );
	    riel            = riel + fe71j;
  	    Jiel(ie71)      = Jiel(ie71) + ge71j;
  	    Hiel(ie71,ie71) = Hiel(ie71,ie71) + He71j;
	 end
	 varargout{1} = riel^2;
	 varargout{2} = 2*Jiel*riel;
	 varargout{3} = 2*(Jiel*Jiel' + riel*Hiel);
      end
   elseif ( iel == (mm+1)*nn+1 )
      switch ( nargout )
      case 1
         riel = -q;
	 for j = 1:mm
	    ie81 = [ 9 j 7 ];
	    ie82 = [ 1 3 ];
	    ie83 = [ mm+j 9 8 ];
	    pe81 = [  1  1 a(j) b(j) c(j) be(j) be1(j) be2(j) ];
	    pe82 = [ bb  0  0    0    0   al(j) al1(j) al2(j) ];
	    pe83 = [ -1 bb  0    0    0   al(j) al1(j) al2(j) ];
            fe81j = hydc20ls( 'exp4prod', x(ie81), pe81 );
            fe82j = hydc20ls( 'poly1prd', x(ie82), pe82 );
            fe83j = hydc20ls( 'poly2prd', x(ie83), pe83 );
	    riel = riel + fe81j + fe82j + fe83j;
	 end
	 varargout{1} = 1.e-10*riel^2;
      case 2
         riel = -q;
	 Jiel = zeros( niel, 1 );
	 for j = 1:mm
	    ie81 = [ 9 j 7 ];
	    ie82 = [ 1 3 ];
	    ie83 = [ mm+j 9 8 ];
	    pe81 = [  1  1 a(j) b(j) c(j) be(j) be1(j) be2(j) ];
	    pe82 = [ bb  0  0    0    0   al(j) al1(j) al2(j) ];
	    pe83 = [ -1 bb  0    0    0   al(j) al1(j) al2(j) ];
            [ fe81j, ge81j ] = hydc20ls( 'exp4prod', x(ie81), pe81 );
            [ fe82j, ge82j ] = hydc20ls( 'poly1prd', x(ie82), pe82 );
            [ fe83j, ge83j ] = hydc20ls( 'poly2prd', x(ie83), pe83 );
	    riel = riel + fe81j + fe82j + fe83j;
	    Jiel(ie81) = Jiel(ie81) +ge81j;
	    Jiel(ie82) = Jiel(ie82) +ge82j;
	    Jiel(ie83) = Jiel(ie83) +ge83j;
	 end
	 varargout{1} = 1.e-10*riel^2;
	 varargout{2} = 2.e-10*Jiel*riel;
      case 3
         riel = -q;
	 Jiel = zeros( niel, 1 );
	 Hiel = zeros( niel, niel );
	 for j = 1:mm
	    ie81 = [ 9 j 7 ];
	    ie82 = [ 1 3 ];
	    ie83 = [ mm+j 9 8 ];
	    pe81 = [  1  1 a(j) b(j) c(j) be(j) be1(j) be2(j) ];
	    pe82 = [ bb  0  0    0    0   al(j) al1(j) al2(j) ];
	    pe83 = [ -1 bb  0    0    0   al(j) al1(j) al2(j) ];
            [ fe81j, ge81j, He81j ] = hydc20ls( 'exp4prod', x(ie81), pe81 );
            [ fe82j, ge82j, He82j ] = hydc20ls( 'poly1prd', x(ie82), pe82 );
            [ fe83j, ge83j, He83j ] = hydc20ls( 'poly2prd', x(ie83), pe83 );
	    riel = riel + fe81j + fe82j + fe83j;
	    Jiel(ie81) = Jiel(ie81) +ge81j;
	    Jiel(ie82) = Jiel(ie82) +ge82j;
	    Jiel(ie83) = Jiel(ie83) +ge83j;
            Hiel(ie81,ie81) = Hiel(ie81,ie81) + He81j;
            Hiel(ie82,ie82) = Hiel(ie82,ie82) + He82j;
            Hiel(ie83,ie83) = Hiel(ie83,ie83) + He83j;
	 end
	 varargout{1} = 1.e-10*riel^2;
	 varargout{2} = 2.e-10*Jiel*riel;
         varargout{3} = 2.e-10*(Jiel*Jiel' + riel*Hiel);
      end
   else
      i = iel - ((mm+1)*nn+1);
      switch ( nargout )
      case 1
         riel = 0;
	 if ( i == kk )
	    riel = - smallHF;
	 elseif ( i == kk+1 )
	    riel = - bigHF;
	 end
         for j = 1:mm
	    ie91 = [    14   mm+j 11 ];
	    ie92 = [   mm+j   13  11 ];
	    ie93 = [    13     j  10 ];
	    ie94 = [ mm+mm+j  14  12 ];
	    pe91 = [  1 1 a(j) b(j) c(j) be(j) be1(j) be2(j) ];
	    if ( i < kk )
  	       pe92 = [  1 bb  0    0    0   al(j) al1(j) al2(j) ];
	       pe94 = [ -1 bb  0    0    0   al(j) al1(j) al2(j) ];
	    elseif ( i == kk )
  	       pe92 = [  1 bb  0    0    0   al(j) al1(j) al2(j) ];
	       pe94 = [ -1 -d  0    0    0   al(j) al1(j) al2(j) ];
	    else
  	       pe92 = [  1 -d  0    0    0   al(j) al1(j) al2(j) ];
	       pe94 = [ -1 -d  0    0    0   al(j) al1(j) al2(j) ];
	    end
	    pe93 = [ -1 1 a(j) b(j) c(j) be(j) be1(j) be2(j) ];
            fe91j = hydc20ls( 'exp4prod', x(ie91), pe91 );
            fe92j = hydc20ls( 'poly2prd', x(ie92), pe92 );
            fe93j = hydc20ls( 'exp4prod', x(ie93), pe93 );
            fe94j = hydc20ls( 'poly2prd', x(ie94), pe94 );
	    riel  = riel + fe91j + fe92j + fe93j + fe94j;
	 end
         varargout{1} = 1.e-10*riel^2;
      case 2
         riel = 0;
         Jiel = zeros( niel, 1 );
         for j = 1:mm 
	    ie91 = [    14   mm+j 11 ];
	    ie92 = [   mm+j   13  11 ];
	    ie93 = [    13     j  10 ];
	    ie94 = [ mm+mm+j  14  12 ];
	    pe91 = [  1 1 a(j) b(j) c(j) be(j) be1(j) be2(j) ];
	    pe92 = [  1 0  0    0    0   al(j) al1(j) al2(j) ];
	    pe93 = [ -1 1 a(j) b(j) c(j) be(j) be1(j) be2(j) ];
	    pe94 = [ -1 0  0    0    0   al(j) al1(j) al2(j) ];
            [ fe91j, ge91j ] = hydc20ls( 'exp4prod', x(ie91), pe91 );
            [ fe92j, ge92j ] = hydc20ls( 'poly2prd', x(ie92), pe92 );
            [ fe93j, ge93j ] = hydc20ls( 'exp4prod', x(ie93), pe93 );
            [ fe94j, ge94j ] = hydc20ls( 'poly2prd', x(ie94), pe94 );
	    riel = riel + fe91j + fe92j + fe93j + fe94j;
	    Jiel(ie91) = Jiel(ie91) + ge91j;
	    Jiel(ie92) = Jiel(ie92) + ge92j;
	    Jiel(ie93) = Jiel(ie93) + ge93j;
	    Jiel(ie94) = Jiel(ie94) + ge94j;
	 end
         varargout{1} = 1.e-10*riel^2;
	 varargout{2} = 2.e-10*Jiel*riel;
      case 3
         riel = 0;
         Jiel = zeros( niel, 1 );
         Hiel = zeros( niel, niel );
         for j = 1:mm
	    ie91 = [  14    mm+j 11 ];
	    ie92 = [ mm+j    13  11 ];
	    ie93 = [  13      j  10 ];
	    ie94 = [ mm+mm+j 14  12 ];
	    pe91 = [  1 1 a(j) b(j) c(j) be(j) be1(j) be2(j) ];
	    pe92 = [  1 0  0    0    0   al(j) al1(j) al2(j) ];
	    pe93 = [ -1 1 a(j) b(j) c(j) be(j) be1(j) be2(j) ];
	    pe94 = [ -1 0  0    0    0   al(j) al1(j) al2(j) ];
            [ fe91j, ge91j, He91j ] = hydc20ls( 'exp4prod', x(ie91), pe91 );
            [ fe92j, ge92j, He92j ] = hydc20ls( 'poly2prd', x(ie92), pe92 );
            [ fe93j, ge93j, He93j ] = hydc20ls( 'exp4prod', x(ie93), pe93 );
            [ fe94j, ge94j, He94j ] = hydc20ls( 'poly2prd', x(ie94), pe94 );
	    riel = riel + fe91j + fe92j + fe93j + fe94j;
	    Jiel(ie91) = Jiel(ie91) + ge91j;
	    Jiel(ie92) = Jiel(ie92) + ge92j;
	    Jiel(ie93) = Jiel(ie93) + ge93j;
	    Jiel(ie94) = Jiel(ie94) + ge94j;
	    Hiel(ie91,ie91) = Hiel(ie91,ie91) + He91j;
	    Hiel(ie92,ie92) = Hiel(ie92,ie92) + He92j;
	    Hiel(ie93,ie93) = Hiel(ie93,ie93) + He93j;
	    Hiel(ie94,ie94) = Hiel(ie94,ie94) + He94j;
	 end
         varargout{1} = 1.e-10*riel^2;
	 varargout{2} = 2.e-10*Jiel*riel;
         varargout{3} = 2.e-10*(Jiel*Jiel' + riel*Hiel);
      end
   end

case '2prod'

   x = varargin{1};
   p = varargin{2};
   varargout{1} = p(1)*x(1)*(x(2)+p(2));
   if ( nargout > 1 )
      varargout{2} = [ p(1)*(x(2)+p(2)); p(1)*x(1) ];
      if ( nargout > 2 )
         varargout{3} = [ 0 p(1); p(1) 0 ];
      end
   end
   
case 'poly1prd'

   x    = varargin{1};
   p    = varargin{2};
   poly =  p(6) + p(7) * x(2) + p(8) * x(2)^2;
   varargout{1} = p(1)*x(1)*poly;
   if ( nargout > 1 )
      dpoly =  p(7)+ 2*p(8)*x(2);
      varargout{2} = [ p(1)*poly; p(1)*x(1)*dpoly ];
      if ( nargout > 2 )
         varargout{3} = [  0  p(1)*dpoly; p(1)*dpoly 2*p(1)*x(1)*p(8) ];
      end
   end

case 'poly2prd'

   x    = varargin{1};
   p    = varargin{2};
   poly = p(6) + p(7)*x(3)+p(8)*x(3)^2;
   varargout{1} = p(1)*x(1)*(p(2)+x(2))*poly;
   if ( nargout > 1 )
      dpoly = p(7)+2*p(8)*x(3);
      varargout{2} = [ p(1)*(p(2)+x(2))*poly;
                       p(1)*x(1)*poly;
                       p(1)*x(1)*(p(2)+x(2))*dpoly ];
      if ( nargout > 2 )
         varargout{3} = [           0               p(1)*poly        p(1)*(p(2)+x(2))*dpoly;
	                       p(1)*poly                0              p(1)*x(1)*dpoly;
			  p(1)*(p(2)+x(2))*dpoly  p(1)*x(1)*dpoly  2*p(1)*x(1)*(p(2)+x(2))*p(8) ];
      end
   end
   
case 'exp2prod'

   x      = varargin{1};
   p      = varargin{2};
   p5x2   = p(5)+x(2);
   exprod = p(1)*p(2)*exp(p(3)+p(4)/p5x2);
   varargout{1} = x(1)*exprod;
   if ( nargout > 1 )
      varargout{2} = [ exprod; -x(1)*exprod*p(4)/(p5x2^2) ];
      if ( nargout > 2 )
         varargout{3} = [          0                  -exprod*p(4)/(p5x2^2);
	                 -exprod*p(4)/(p5x2^2)  varargout{1}*((p(4)/(p5x2)^2)^2 + 2*p(4)/(p5x2^3)) ];
      end
   end
   
case 'exp3prod'

   x      = varargin{1};
   p      = varargin{2};
   p5x3   = p(5)+x(3);
   exprod = p(1)*p(2)*exp(p(3)+p(4)/p5x3);
   f      = x(1)*x(2)*exprod;
   varargout{1} = f;
   if ( nargout > 1 )
      trm = -p(4)/(p5x3^2);
      varargout{2} = [ x(2)*exprod; x(1)*exprod; f*trm ];
      if ( nargout > 2 )
         varargout{3} = [   0                    exprod           x(2)*exprod*trm;
	                  exprod                  0               x(1)*exprod*trm;
			  x(2)*exprod*trm   x(1)*exprod*trm   f*(trm^2+2*p(4)/(p5x3^3)) ];
      end
   end

case 'exp4prod'

   x      = varargin{1};
   p      = varargin{2};
   p5x3   = p(5)+x(3);
   exprod = p(1)*p(2)*exp(p(3)+p(4)/p5x3);
   f      = x(1)*x(2)*exprod;
   poly   = p(6)+p(7)*x(3)+p(8)*x(3)^2;
   varargout{1} = f*poly;
   if ( nargout > 1 )
      dpoly = p(7)+2*p(8)*x(3);
      trm   = dpoly - poly * p(4) / (p5x3)^2;
      varargout{2} = [ x(2)*exprod*poly;
                       x(1)*exprod*poly;
		       f * trm ];
      if ( nargout > 2 )
         H33   = f*(-(p(4)/(p5x3^2))*trm   ...
	         +2*p(8)-dpoly*p(4)/(p5x3^2) ...
		 +2*poly*p(4)/(p5x3^3));
         varargout{3} = [   0               exprod*poly      x(2)*exprod*trm;
	                 exprod*poly             0           x(1)*exprod*trm;
                         x(2)*exprod*trm  x(1)*exprod*trm  H33             ];
      end
   end

end

return

end
