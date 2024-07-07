%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = schmvett( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The Schmidt and Vetters problem.
%
%   Source:
%      J.W. Schmidt and K. Vetters,
%      "Albeitungsfreie Verfahren fur Nichtlineare Optimierungsproblem",
%      Numerische Mathematik 15:263-282, 1970.
%   Also problem 14 (p. 90) in 
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%   and problem 35 in
%      Ph.L. Toint,
%      "Test problems for partially separable optimization and results
%      for the routine PSPMIN",
%      Report 83/4, Department of Mathematics, FUNDP (Namur, B), 1983.
%
%   Ph. Toint 22 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'schmvett';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 0.5; 0.5; 0.5];             % x0
   varargout{2} = -3;                           % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'OUR2-AN-3-0';                % class

case 'cpsstr'

   n  = varargin{1};
   eldom = cell( n-2, 1 );
   for iel = 1:n-2
      eldom{ iel } = [ iel iel+1 iel+2 ];
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
   varargout = opm_eval_cpsf( pname, 'elobjf', x, cpsstr.eldom, nargout, {} );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   x     = varargin{2};
   x1mx2 = (x(1)-x(2));
   den   = 1+x1mx2^2;
   A     = -1/den;
   sarg  =  0.5*( pi*x(2)+x(3));
   B     = -sin(sarg);
   C     = (x(1)+x(3))/x(2)-2;
   D     = -exp(-C^2);
   varargout{1} = A+B+D;
   if ( nargout > 1 )
      gA = (2*x1mx2/den^2)*[ 1; -1; 0 ];
      gB = -cos(sarg)*[ 0; 0.5*pi; 0.5];
      gC = [ 1/x(2); -(x(1)+x(3))/x(2)^2; 1/x(2) ];
      gD = -2*D*C*gC;
      varargout{2} = gA+gB+gD;
      if ( nargout > 2 )
         HA = 2*[ 1/den^2*(1-4*x1mx2^2/den), -1/den^2*(1-4*x1mx2^2/den), 0;
	         -1/den^2*(1-4*x1mx2^2/den),  1/den^2*(1-4*x1mx2^2/den), 0;
		            0              ,            0              , 0];
         HB = -B*[ 0; 0.5*pi; 0.5]*[ 0; 0.5*pi; 0.5].';
	 HC = (1/x(2)^2) * [    0,         -1        ,  0;
	                       -1, 2*(x(1)+x(3))/x(2), -1;
	                        0,         -1        ,  0 ];
         HD = 2*D*((2*C^2-1)*gC*gC.'-C*HC);
         varargout{3} = HA+HB+HD;
      end
   end

end

return

end