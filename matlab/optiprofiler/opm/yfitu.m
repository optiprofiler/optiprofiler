%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = yfitu( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   A nonlinear least-squares problem.  This problem arises in measuring
%   angles and distances to a vibrating beam using a laser-Doppler
%   velocimeter.

%   Source:
%   an exercize for L. Watson course on LANCELOT in the Spring 1993.

%   Ph. Toint, 25 VII 2021, from the CUTEst YFITU.SIF
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'yfitu';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 0.6; -0.6; 20.0 ];          % x0
   varargout{2} = 0.0;                          % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-MN-3-0';                % class

case 'cpsstr'

   y = [ 21.158931;   17.591719;   14.046854;  10.519732;
          7.0058392;   3.5007293;   0.000000;  -3.5007293;
	 -7.0058392; -10.519732;  -14.046854; -17.591719;
	-21.158931;  -24.753206;  -28.379405; -32.042552; -35.747869 ];
   cpsstr.name  = pname;
   cpsstr.eldom = { [ 1 2 3 ] };
   cpsstr.param = { y };
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

   x    = varargin{2};
   y    = varargin{3};
   varargout{1} = 0;
   varargout{2} = zeros(3,1);
   varargout{3} = zeros(3,3);
   for i = 1:17
     f1    = (i-1)/16;
     ttan  = tan( x(1)*(1-f1) + x(2)*f1 );
     tsec2 = 1/cos( x(1)*(1-f1) + x(2)*f1 )^2;
     ri    =  x(3) * ttan - y(i);
     varargout{1} =  varargout{1} + ri^2;
     if ( nargout > 1 )
        Ji = [ x(3)*(1-f1)*tsec2; x(3)*f1*tsec2; ttan ];
	varargout{2} = varargout{2} + 2 * Ji * ri;
	if ( nargout > 2 )
           Hi = [ 2*x(3)*((1.0-f1)^2)*tsec2*ttan   2*x(3)*(1.0-f1)*f1*tsec2*ttan   (1.0-f1)*tsec2;
	          2*x(3)*(1.0-f1)*f1*tsec2*ttan    2*x(3)*(f1^2)*tsec2*ttan        f1*tsec2;
		  (1.0-f1)*tsec2                   f1*tsec2                        0              ];
           varargout{3} = varargout{3} + 2 * ( Ji*Ji' + ri*Hi );
	end
     end
   end
end

return

end
   
