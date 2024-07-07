%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = spmsqrt( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The tridiagonal matrix square root problem by Nocedal and Liu.
%   seen as a nonlinear least-squares problem.

%   Source:  problem 151 (p. 92-93) in
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

pname   = 'spmsqrt';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      m = (n+2)/3;
      if ( n < 1 || round( m ) ~= m )
         disp( [ ' ERROR in spmsqrt: n = ', int2str(n),' but should satisfy n >= 4, n = 3m-2!' ] )
      end
   else
      n = 10;
      m = 4;
   end
   varargout{1} = 0.2*sin( ([1:n].^2)' );       % x0
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AY-V-0';                % class

case 'cpsstr'

   %  Define one least-square residual for each entry of the nel entries of the
   %  column-ordered A.
   %  iv is the integer vector such that each product of index iel is
   %       x(eldom{iel}(iv{iel}(1:l2))).'*x(eldom{iel}(iv{iel}(l2+1:l)))
   %  where l = length(iv{iel}} and l2 = l/2.  This trick is avoids repetitions
   %  in eldom{iel}.

   n     = varargin{1};
   m     = (n+2)/3;
   nel   = 5*m-6;
   rows  = cell( m, 1 );        % the variables occuring in each row
   cols  = cell( m, 1 );        % the variables occuring in each column
   for j = 1:m
      if ( j == 1 )
         cols{ 1 } = [ 1 2 ];
	 rows{ 1 } = [ 1 3 ];
      elseif ( j == m )
         cols{ m } = [ n-1 n ];
	 rows{ m } = [ n-2 n ];
      else
         cols{ j } = [ 3*(j-1):3*j-1 ];
	 rows{ j } = [ 3*(j-1)-1, 3*(j-1)+1, 3*(j-1)+3 ];
      end
   end
   eldom = cell( nel, 1 );
   iv    = cell( nel, 1 );
   iel   = 0;
   for j = 1:m
       for i = max(1,j-2):min(j+2,m)
         iel = iel + 1;
	 switch ( iel )
	 case 1
	    eldom{ 1 } = [ 1 2 3 ];
	    iv   { 1 } = [ 1 3 1 2 ];
	 case 2
	    eldom{ 2 } = [ 1 2 4 ];
	    iv   { 2 } = [ 2 3 1 2 ];
	 case 3
	    eldom{ 3 } = [ 2 5 ];
	    iv   { 3 } = [ 2 1 ];
	 case 4
	    eldom{ 4 } = [ 1 3 4 ];
	    iv   { 4 } = [ 1 2 2 3 ];
	 case 8  
	    eldom{ 8 } = [ 3 6 ];
	    iv   { 8 } = [ 1 2 ];
	 case nel-7
	    eldom{ iel } = [ n-2 n-5 ];
	    iv   { iel } = [ 1 2 ];
	 case nel-3
	    eldom{ iel } = [ n-2 n n-3 ];
	    iv   { iel } = [ 1 2 3 1 ];
	 case nel-2
	    eldom{ iel } = [ n-4 n-1 ];
	    iv   { iel } = [ 1 2 ];
	 case nel-1
	    eldom{ iel } = [ n-3 n-1 n ];
	    iv   { iel } = [ 1 2 2 3 ];
	 case nel
	    eldom{ nel } = [ n-2 n n-1 ];
	    iv   { nel } = [ 1 2 3 2 ];
	 otherwise
    	    if ( i == j-2 )
	        eldom{ iel } = [ rows{i}(3) cols{j}(1) ];
		iv   { iel } = [ 1 2 ];
	    elseif ( i == j-1 )
	        eldom{ iel } = [ rows{i}(2:3) cols{j}(2) ];
		iv   { iel } = [ 1 2 2 3 ];
	    elseif ( i == j )
                eldom{ iel } = [ rows{i} cols{j}(1) cols{j}(3) ];
		iv   { iel } = [ 1 2 3 4 2 5 ];
	    elseif ( i == j+1 )
	        eldom{ iel } = [ rows{i}(1:2) cols{j}(2) ];
		iv   { iel } = [ 1 2 3 1 ];
	    else
	        eldom{ iel } = [ rows{i}(1) cols{j}(3) ];
		iv   { iel } = [ 1 2 ];
            end
	 end
      end
   end
   
%  Define B (columnwise)
   
   B     = sparse( m, m );
   b     = sin( ([1:3*m-2].^2).' );
   for j = 1:m
      if ( j == 1 )
	 B( [ 1,2 ], 1 )     = b(1:2);
	 ib = 2;
      elseif ( j == m )
	 B( [ m-1,m ], m )   = b( ib+1:ib+2 );
      else
	 B( [j-1,j,j+1], j ) = b(ib+1:ib+3);
	 ib = ib + 3;
      end
   end

%  Define A

   A = B*B;
   a = [];
   for j = 1:m
      a = [ a; A(max(1,j-2):min(j+2,m),j) ];
   end

   cpsstr.name  = pname;
   cpsstr.eldom = eldom;
   cpsstr.param = { a, iv };
   varargout{1} = cpsstr;

case 'objf'   % varargout = [ f, g, H ]

   x = varargin{1};
   if ( nargin > 2 && isstruct( varargin{2} ) && isfield( varargin{2}, 'eldom' ) )
      cpsstr = varargin{2};
   else
      cpsstr = problem( 'cpsstr', length( x ) );
   end
   varargout = opm_eval_cpsf( pname, 'elobjf', x, cpsstr.eldom, nargout, cpsstr.param{:} );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   iel  = varargin{1};
   x    = varargin{2};
   lx   = length( x );
   iv   = varargin{4};
   liv  = length( iv{iel} );
   liv2 = liv / 2;
   ievr = iv{iel}(1:liv2);
   ievc = iv{iel}(liv2+1:liv);
   riel = varargin{3}(iel) - x(ievr).'*x(ievc);
   varargout{1} = riel^2;
   if ( nargout > 1 )
      Jiel  = zeros( lx, 1 );
      for k = 1:liv2
         Jiel(ievr(k))= Jiel(ievr(k))-x(ievc(k));
	 Jiel(ievc(k))= Jiel(ievc(k))-x(ievr(k));
      end
      varargout{2} = 2*Jiel*riel;
      if ( nargout > 2 )
         Hiel = sparse( lx, lx );
         for k = 1:liv2
	    Hiel(ievc(k),ievr(k)) = Hiel(ievc(k),ievr(k))-1;
	    Hiel(ievr(k),ievc(k)) = Hiel(ievr(k),ievc(k))-1;
	 end
         varargout{3} = 2* (Jiel*Jiel.'+riel*Hiel);
     end
   end

end

return

end