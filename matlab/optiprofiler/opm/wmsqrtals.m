%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = wmsqrtals( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The dense matrix square root problem by Nocedal and Liu.
%   (Case 0) seen as a nonlinear least-squares problem.  This problem
%   results from an error in specifying the msqrtals problem.

%   Source for msqrtals: problem 201 (p. 93) in
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

pname   = 'wmsqrtals';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n     = varargin{1};
      sqrtn = sqrt( n );
      if ( n < 1 || round( sqrtn ) ~= sqrtn )
         disp( [ ' ERROR in wmsqrtals: n = ', int2str(n),' but should satisfy n >= 4, n = m^2!' ] )
      end
   else
      n = 16;
   end
   varargout{1} = 0.2*sin( ([1:n].^2)' );       % x0
   varargout{2} = 'unknown';                    % fstar
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
         eldom{ iel } = [ iel, setdiff(irow,iel), setdiff(icol,iel) ];
	 iel          = iel + 1;
      end
   end
   b = sin( ([1:n].^2).' );
   B = reshape( b, sqrtn, sqrtn );
   A = B*B;
   cpsstr.name  = pname;
   cpsstr.eldom = eldom;
   cpsstr.param = { A(:) };
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

   iel   = varargin{1};
   x     = varargin{2};
   lx    = length( x );
   sqrtn = (lx+1)/2;
   riel  = varargin{3}(iel) - x(1)^2 - x(2:sqrtn).'*x(sqrtn+1:lx);
   varargout{1} = riel^2;
   if ( nargout > 1 )
      Jiel             =  zeros( lx, 1 );
      Jiel(1)          = -2*x(1);
      Jiel(2:sqrtn)    = -x(sqrtn+1:lx);
      Jiel(sqrtn+1:lx) = -x(2:sqrtn);
      varargout{2} = 2*Jiel*riel;
      if ( nargout > 2 )
         Hiel      = sparse( lx, lx );
	 Hiel(1,1) = -2;
         for k = 1:sqrtn-1
	     Hiel(1+k,sqrtn+k) = -1;
	     Hiel(sqrtn+k,1+k) = -1;
	 end
         varargout{3} = 2* (Jiel*Jiel.'+riel*Hiel);
     end
   end
end

return

end