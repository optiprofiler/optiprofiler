%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = tridia( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Shanno's TRIDIA quadratic tridiagonal problem
%
%   Source: problem 8 in
%      Ph.L. Toint,
%      "Test problems for partially separable optimization and results
%      for the routine PSPMIN",
%      Report 83/4, Department of Mathematics, FUNDP (Namur, B), 1983.
%   Also problem 40 (p. 96) in
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

pname   = 'tridia';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 2 )
         disp( [ ' ERROR in tridia: n = ', int2str(n),' should be > 1!' ] )
      end
   else
      n = 10;
   end
   varargout{1} = ones( n, 1);                  % x0
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'QUR2-AY-V-0';                % class

case 'cpsstr'

   n          = varargin{1};
   eldom      = cell( n, 1 );
   eldom{ 1 } = [ 1 ];
   for iel = 2:n
      eldom{ iel } = [ iel-1 iel ];
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

   iel   = varargin{1};
   x     = varargin{2};
   if ( iel == 1 )
      riel = x(1) - 1;
      varargout{1} = riel^2;
      if ( nargout > 1 )
         varargout{2} =  2 * riel;
	 if ( nargout > 2 )
	    varargout{3} = 2;
	 end
      end
   else
      riel  = 2 * x(2) - x(1);
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = [ -1; 2 ];
         varargout{2} = 2 * Jiel * riel;
         if ( nargout > 2 )
            varargout{3} = 2 * Jiel * Jiel.';
         end
      end
   end

return

end