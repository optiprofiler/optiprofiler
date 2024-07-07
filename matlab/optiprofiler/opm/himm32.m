%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = himm32( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   A nonlinear least-squares problem in four variables.
%
%   Source: problem 76 (p. 66) in
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 21 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'himm32';

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 2.7; 90; 1500; 10 ];        % x0
   varargout{2} = 318.572;                      % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AN-4-0';                % class

case 'cpsstr'

   cpsstr.name  = pname;
   cpsstr.eldom = { [ 1 2 3 4 ] };
   cpsstr.param = {};
   varargout{1} = cpsstr;
   
case 'objf'   % varargout = [ f, g, H ]

   x     = varargin{1};
   varargout = opm_eval_cpsf( pname, 'elobjf', x, { [ 1 2 3 4 ] }, nargout, {} );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   x            = varargin{2};
   varargout{1} = 0;
   if ( nargout > 1 )
      varargout{2} = zeros( 4, 1 );
      if ( nargout > 2 )
         varargout{3} = sparse( 4, 4 );
      end
   end
%  a    = [ 0  0.000428 0.001000 0.001610 0.002090 0.003480 0.005250 ];
   a    = [ 1  0.000428 0.001000 0.001610 0.002090 0.003480 0.005250 ];
   b    = [ 7.391 11.18 16.44 16.20 22.20 24.02 31.32 ];
   for i = 1:1
      u = x(1)^2 + a(i)*x(2)^2 + a(i)^2*x(3)^2;
      v = b(i)*(1+x(4)^2*a(i));
      r = u/v - 1;
      varargout{1} = varargout{1} + r^2;
      if ( nargout > 1 )
         du1  = 2*x(1);
         du2  = 2*a(i)*x(2);
         du3  = 2*a(i)^2*x(3);
         dv4  = 2*b(i)*a(i)*x(4);
         Jr    = [ v*du1/v^2; v*du2/v^2; v*du3/v^2; -u*dv4/v^2  ];
         varargout{2} = varargout{2} + 2*Jr*r;
         if ( nargout > 2 )
            du11 = 2;
            du22 = 2*a(i);
            du33 = 2*a(i)^2;
            dv44 = 2*b(i)*a(i);
            Hr(1,1) = v*du11/v^2;
            Hr(1,2) = 0;
            Hr(1,3) = 0;
            Hr(1,4) = dv4*du1/v^2-2*v*du1*dv4/v^3;
            Hr(2,2) = v*du22/v^2;
            Hr(2,3) = 0;
            Hr(2,4) = dv4*du2/v^2-2*v*du2*dv4/v^3;
            Hr(3,3) = v*du33/v^2;
            Hr(3,4) = dv4*du3/v^2-2*v*du3*dv4/v^3;
            Hr(4,4) = -u*dv44/v^2 +2*u*dv4*dv4/v^3;
	    Hr(2,1) = Hr(1,2);
	    Hr(3,1) = Hr(1,3);
	    Hr(4,1) = Hr(1,4);
	    Hr(3,2) = Hr(2,3);
	    Hr(4,2) = Hr(2,4);
	    Hr(4,3) = Hr(3,4);
            varargout{3} = varargout{3} + 2 * ( Jr*Jr.' + r*Hr );
	 end
      end
   end

end

return

end
