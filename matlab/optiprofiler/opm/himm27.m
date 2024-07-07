%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = himm27( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   A nonlinear problem in two variables with minimizers at infinity.
%   It is similar to a nonlinear-least-squares problem, except that
%   the objective function is the product of the squared residuals instead of
%   their sum.
%
%   Source: problem  77 (p. 62) in
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 21 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ -1.2; 1 ];                  % x0
   varargout{2} = 0        ;                    % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'OUR2-AN-2-0';                % class

case 'cpsstr'

   cpsstr.name  ='himm27';
   cpsstr.eldom = { [ 1 2 ] };
   cpsstr.param = {};
   varargout{1} = cpsstr;
   
case 'objf'   % varargout = [ f, g, H ]

   x     = varargin{1};
   varargout = opm_eval_cpsf( 'himm27', 'elobjf', x, { [ 1 2 ] }, nargout, {} );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   x    = varargin{2};
   r1   = x(1)*x(2);
   r2   = 1-x(1);
   r3   = 1 - x(2) - x(1)*(1-x(1))^5;
   varargout{1} = r1^2 * r2^2 * r3^2;
   if ( nargout > 1 )
      gr1  = [ x(2); x(1) ];
      gr2  = [ -1; 0 ];
      gr3  = [ (1-x(1))^4*(6*x(1)-1); -1 ];
      g    =  gr1*r2*r3 + r1*gr2*r3 + r1*r2*gr3;  
      varargout{2} = 2*r1*r2*r3*g;
      if ( nargout > 2 )
         Hr1 = [ 0, 1; 1, 0 ];
	 Hr2 = [ 0, 0; 0, 0 ];
	 Hr3 = [ -4*(1-x(1))^3*(6*x(1)-1)+6*(1-x(1))^4, 0;
	                     0,                         0];
	 H   = Hr1*r2*r3    + gr1*gr2.'*r3  + gr1*r2*gr3.'+...
	       gr2*r3*gr1.' + r1*Hr2*r3     + r1*gr2*gr3.'+ ...
	       r2*gr3*gr1.' + r1*gr3*gr2.'  + r1*r2*Hr3   ;
         varargout{3} = 2*g*g.'+2*r1*r2*r3*H;
      end
   end
end

return

end
