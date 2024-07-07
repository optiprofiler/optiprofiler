%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = jensmp( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The Jennrich and Sampson problem
%
%   Source:  Problem 6 in
%      J.J. More', B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%   Also problem 23 (p. 69) in 
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 21 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname = 'jensmp';

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 0.3; 0.4];                  % x0
   varargout{2} = 124.362;                      % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AN-2-0';                % class

case 'cpsstr'

   cpsstr.name  = pname;
   cpsstr.eldom = { [ 1 2 ] };
   cpsstr.param = {};
   varargout{1} = cpsstr;
   
case 'objf'   % varargout = [ f, g, H ]

   varargout = opm_eval_cpsf( pname, 'elobjf', varargin{1}, { [ 1 2 ] }, nargout, {} );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   x    = varargin{2};
   varargout{1} = 0;
   if ( nargout > 1 )
      varargout{2} = zeros(2,1);
      if ( nargout > 2 )
         varargout{3} = zeros( 2, 2 );
      end
   end
   for iel = 1:10
      e1   = exp(iel*x(1));
      e2   = exp(iel*x(2));
      riel = 2 + 2*iel - e1 - e2;
      varargout{1} = varargout{1} + riel^2;
      if ( nargout > 1 )
         Jiel = iel * [ -e1; -e2 ];
         varargout{2} = varargout{2} + 2* Jiel * riel;
         if ( nargout > 2 )
            Hiel = iel^2 *[ -e1, 0;
	              0, -e2 ];
	    varargout{3} = varargout{3} + 2 * ( Jiel*Jiel.' + riel*Hiel );
         end
      end
   end
end

return

end