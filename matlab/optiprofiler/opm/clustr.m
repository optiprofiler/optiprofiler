%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = clustr( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The Clustr nonlinear least-squares problem in 2 variables.
%
%   Source: problem 207 in 
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 18 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'clustr';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 0; 0];                      % x0
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AN-2-0';                % class

case 'cpsstr'

   eldom{ 1 }   = [ 1 2 ];
   eldom{ 2 }   = [ 1 2 ];
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

   iel  = varargin{1};
   x    = varargin{2};
   switch ( iel )
   case 1
      t1   = x(1)-x(2)^2;
      t2   = x(1) - sin( x(2) );
      riel = t1 * t2;
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jt1  = [ 1;  -2*x(2)   ];
	 Jt2  = [ 1; -cos(x(2)) ];
	 Jiel = t2*Jt1 + t1*Jt2;
	 varargout{2} = 2 * Jiel * riel;
	 if ( nargout > 2 )
            Ht1 = [ 0,    0;
	            0,   -2      ];
            Ht2 = [ 0,    0;
	            0, sin(x(2)) ];
            Hiel = Jt2*Jt1.'+t2*Ht1 + Jt1*Jt2.'+t1*Ht2;
	    varargout{3} = 2 * ( Jiel*Jiel.' + riel*Hiel );
	 end
      end
   case 2
      t1 = cos(x(2)) - x(1);
      t2 = x(2) - cos(x(1));
      riel = t1 * t2;
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jt1 = [ -1; -sin(x(2))];
         Jt2 = [ sin(x(1)); 1  ];
	 Jiel = t2*Jt1 + t1*Jt2;
	 varargout{2} = 2 * Jiel * riel;
          if ( nargout > 2 )
             Ht1 = [  0,     0;
	              0, -cos(x(2)) ];
	     Ht2 = [ cos(x(1)), 0;
	                0,      0 ];
             Hiel = Jt1*Jt2.' + t2*Ht1 + Jt2*Jt1.' + t1*Ht2;
	     varargout{3} = 2 * ( Jiel*Jiel.' + riel*Hiel );
	  end
       end
   end
end

return

end