%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = kowosb( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The Kowalik and Osborne data-fitting problem.
%   A problem arising in the analysis of kinetic data for an enzyme
%   reaction.
%
%   Source:  Problem 15 in
%      J.J. More', B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%   Also problem 31 (p. 70) in 
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 21 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname = 'kowosb';

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 0.25; 0.39; 0415; 0.39 ];   % x0
   varargout{2} = 0.00307505;                   % fstar
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

   varargout = opm_eval_cpsf( pname, 'elobjf', varargin{1}, { [ 1 2 3 4 ] }, nargout, {} );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   x    = varargin{2};
   varargout{1} = 0;
   if ( nargout > 1 )
      varargout{2} = zeros( 4, 1 );
      if ( nargout > 2 )
         varargout{3} = zeros( 4, 4 );
      end
   end
   y = [ 0.1957 0.1947 0.1735 0.1600 0.0844 0.0627 0.0456 0.0342  0.0323 0.0235 0.0246 ];
   u = [ 4.0 2.0 1.0 0.5 0.25 0.167 0.125 0.1 0.0833 0.0714 0.0624 ];
   for iel = 1:1
      ui   = u(iel);
      b1   = ui^2+ui*x(2);
      b2   = ui^2+ui*x(3)+x(4);
      riel = x(1)*b1/b2-y(iel);
      varargout{1} = varargout{1} + riel^2;
      if ( nargout > 1 )
         ux1  = ui*x(1);
	 b2sq = b2^2;
	 t1   = b1/b2sq;
         Jiel = [ b1/b2; ux1/b2; -ux1*t1; -x(1)*t1 ];
         varargout{2} = varargout{2} + 2* Jiel * riel;
         if ( nargout > 2 )
	    ub1  = ui*b1;
	    t2   = 2/b2^3;
            Hiel = [     0    ,      ui/b2    , -ub1/b2sq ,    -t1    ;
	               ui/b2  ,        0    , -ux1*ui/b2sq, -ux1/b2sq ;
		     -ub1/b2sq, -ux1*ui/b2sq,  t2*ux1*ub1 , t2*ux1*b1 ;
		        -t1   ,    -ux1/b2sq,  t2*ux1*b1  , t2*x(1)*b1 ];
	    varargout{3} = varargout{3} + 2 * ( Jiel*Jiel.' + riel*Hiel );
         end
      end
   end
end

return

end