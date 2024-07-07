%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = s308( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Source: problem 308 in
%   K. Schittkowski,
%   " More Test Problems for Nonlinear Programming Codes",
%   Springer Verlag, Berlin, 1987.
%
%   SIF input: Ph. Toint, April 1991.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 's308';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 3; 0.1];                     % x0
   varargout{2} = 0.773199;                       % fstar
   varargout{3} = '';                             % xtype
   varargout{4} = [];                             % xlower
   varargout{5} = [];                             % xupper
   varargout{6} = [];                             % clower
   varargout{7} = [];                             % cupper
   varargout{8} = 'SUR2-AN-2-0';                  % class

case 'cpsstr'

   eldom{1}     = [ 1 2 ];
   eldom{2}     = [ 1 ];
   eldom{3}     = [ 2 ];
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

   iel   = varargin{1};
   x     = varargin{2};
   switch( iel )
   case 1
      riel = x(1)^2 + x(1)*x(2) + x(2)^2;
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = [ 2*x(1)+x(2); x(1)+2*x(2) ];
	 varargout{2} = 2*Jiel*riel;
	 if ( nargout > 2 )
	    Hiel = [ 2   1;   1   2 ];
	    varargout{3} = 2 * ( Jiel*Jiel' + riel*Hiel );
	 end
      end	 
   case 2
      riel = sin( x(1) );
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = [ cos( x(1) ) ];
	 varargout{2} = 2*Jiel*riel;
	 if ( nargout > 2 )
	    Hiel = [ -riel ];
	    varargout{3} = 2 * ( Jiel*Jiel' + riel*Hiel );
	 end
      end	 
   case 3
      riel = cos( x(1) );
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = [ -sin( x(1) ) ];
	 varargout{2} = 2*Jiel*riel;
	 if ( nargout > 2 )
	    Hiel = [ -riel ];
	    varargout{3} = 2 * ( Jiel*Jiel' + riel*Hiel );
	 end
      end	 
   end
end

return

end
