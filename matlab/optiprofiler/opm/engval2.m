%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = engval2( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The ENGVAL2 problem.
%
%   Source: problem 15 (p. 53) in
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint,  19 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'engval2';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 1; 2; 0 ];                  % x0
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AY-V-0';                % class

case 'cpsstr'

   for iel = 1:5
      eldom{ iel } = [ 1 2 3 ];
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

   iel = varargin{1};
   x   = varargin{2};
   switch ( iel )
   case 1
      riel = x(1)^2 + x(2)^2 + x(3)^2 - 1;
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = 2 * x;
	 varargout{2} = 2 * Jiel * riel;
	 if ( nargout > 2 )
            varargout{3} = 2 * ( Jiel*Jiel.' + 2*riel*eye(3));
	 end
      end
   case 2
      riel = x(1)^2 + x(2)^2 + (x(3)-2)^2 - 1;
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = 2 * [ x(1); x(2); x(3)-2 ];
	 varargout{2} = 2 * Jiel * riel;
	 if ( nargout > 2 )
            varargout{3} = 2 * ( Jiel*Jiel.' + 2*riel*eye(3));
	 end
      end
   case 3
      riel = x(1) + x(2) + x(3) - 1;
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = [ 1; 1; 1 ];
	 varargout{2} = 2 * Jiel * riel;
	 if ( nargout > 2 )
            varargout{3} = 2 * Jiel*Jiel.';
	 end
      end
   case 4
      riel = x(1) + x(2) - x(3) - 1;
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = [ 1; 1; -1 ];
	 varargout{2} = 2 * Jiel * riel;
	 if ( nargout > 2 )
            varargout{3} = 2 * Jiel*Jiel.';
	 end
      end
   case 5
      riel = x(1)^3 + 3*x(2)^2 +(5*x(3)-x(1) +1)^2 -36;
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = [ 3*x(1)^2-2*(5*x(3)-x(1) +1); 6*x(2); 10*(5*x(3)-x(1) +1)];
	 varargout{2} = 2 * Jiel * riel;
         if ( nargout > 2 )
            Hiel = [ 6*x(1)+2, 0, -10;
	               0     , 6,   0;
		     -10     , 0,  50 ];
           varargout{3} = 2 * ( Jiel*Jiel.' + riel*Hiel );
	 end
      end
   end
   
return

end