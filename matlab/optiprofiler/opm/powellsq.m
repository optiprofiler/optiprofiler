%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = powellsq( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Source:
%      M.J.D. Powell,
%      " A hybrid method for nonlinear equations",
%      In P. Rabinowitz(ed.) "Numerical Methods for Nonlinear Algebraic
%      Equations", Gordon and Breach, 1970.
%   Also problem 217 (p.84.)
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 24 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname = 'powellsq';

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 3; 1];                      % x0
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'OUR2-AN-2-0';                % class

case 'cpsstr'

   cpsstr.name  = pname;
   cpsstr.eldom = { [ 1 2 ] };
   cpsstr.param = {};
   varargout{1} = cpsstr;
   
case 'objf'   % varargout = [ f, g, H ]

   varargout = opm_eval_cpsf( 'powellsq', 'elobjf', varargin{1}, { [ 1 2 ] }, nargout, {} );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   x  = varargin{2};
   r1 = x(1);
   r2 = 10*x(1)/(x(1)+0.1) + 2*x(2)^2;
   varargout{1} = r1^2 + r2^2;
   if ( nargout > 1 )
      J1 = [ 1; 0 ];
      J2 = [ 10/(x(1)+0.1)-10*x(1)/(x(1)+0.1)^2; 4*x(2) ];
      varargout{2} = 2 * ( J1 * r1 + J2 * r2 );
      if ( nargout > 2 )
         H2 = [ -20/(x(1)+0.1)^2+20*x(1)/(x(1)+0.1)^3, 0;
	                       0                     , 4 ];
         varargout{3} = 2 * ( J1 * J1.' +J2 * J2.' + r2 * H2 );
      end
   end
end

return

end