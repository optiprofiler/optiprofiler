%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = msqrtbls( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The dense matrix square root problem by Nocedal and Liu.
%   (Case 1) seen as a nonlinear least-squares problem.

%   Source:  problem 204 (p. 93) in
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   If the dimension is unspecified, the default n = 16 is chosen.
%
%   Ph. Toint 22 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'msqrtbls';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n     = varargin{1};
      sqrtn = sqrt( n );
      if ( n < 1 || round( sqrtn ) ~= sqrtn )
         disp( [ ' ERROR in msqrtbls: n = ', int2str(n),' but should satisfy n >= 4, n = m^2' ] )
      end
   else
      n = 16;
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

   n     = varargin{1};
   sqrtn = sqrt(n);
   eldom = cell( n, 1 );
   iel   = 1;
   for j = 1:sqrtn
      icol  = [1:sqrtn]+(j-1)*sqrtn;
      for i = 1:sqrtn
         irow         = [0:sqrtn-1]*sqrtn+i;
         eldom{ iel } = [ irow, setdiff(icol,iel) ];
	 iel          = iel + 1;
      end
   end
   b            = sin( ([1:n].^2).' );
   b(2*sqrtn+1) = 0;                                    % Defines Case 1
   B            = reshape( b, sqrtn, sqrtn )';
   A            = B*B;
   cpsstr.name  = pname;
   cpsstr.eldom = eldom;
   cpsstr.param = { A(:) };
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

   %  In what follows:
   %  iel  is the vector index of the current (i,j) pair, that is its index in the vector version of the matrix
   %  d    is the matrix dimension
   %  irow is the ordered list of vector indeces of entries in row i
   %  icol is the ordered list of vector indeces of entries in column j
   %  pir  is the position of iel in irow
   %  pic  is the position of iel in icol
   %  ievr is the position if the numerical values of vector indeces irow in the argument x
   %  ievc is the position if the numerical values of vector indeces icol in the argument x
   
   iel = varargin{1};
   x   = varargin{2};
   lx  = length( x );
   d   = ( lx + 1 ) / 2;            %  the matrix dimension
   j   = floor( iel / d ) + 1;      %  find the i and j corresponding to the vector index iel
   i   = mod( iel, d );
   if ( i == 0 )
      i = d;
      j = j-1;
   end
   irow = [0:d-1]*d+i;              %  build irow
   pir  = find( irow == iel );      %  find iel in it
   icol = [1:d]+(j-1)*d;            %  build icol
   pic  = find( icol == iel );      %  find iel in it
   ievr = [ 1:d ];                  %  find the x values corresponding to row entries
   if ( pic == 1 )                  %  reinset the index corresponding to the missing iel ...
      ievc = [ pir, d+1:lx ];       %  ... in x(d+1;lx) in order to find the x values for ...
   elseif ( pic == d )              %  ... the column
      ievc = [ d+1:lx, pir ];
   else
      ievc = [ d+1:d+pic-1, pir, d+pic:lx ];
   end
   riel   = varargin{3}(iel) - x(ievr).'*x(ievc);
   varargout{1} = riel^2;
   if ( nargout > 1 )
      Jiel =  zeros( lx, 1 );
      for k = 1:d
         Jiel(ievr(k))= Jiel(ievr(k))-x(ievc(k));
	 Jiel(ievc(k))= Jiel(ievc(k))-x(ievr(k));
      end
      varargout{2} = 2*Jiel*riel;
      if ( nargout > 2 )
         Hiel = sparse( lx, lx );
         for k = 1:d
	    Hiel(ievc(k),ievr(k)) = Hiel(ievc(k),ievr(k))-1;
	    Hiel(ievr(k),ievc(k)) = Hiel(ievr(k),ievc(k))-1;
	 end
         varargout{3} = 2* (Jiel*Jiel.'+riel*Hiel);
     end
   end
end

return

end