%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = zangwil3( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The Zangwill quadratic problem in 3 variables.
%
%   Source: problem 13 (p. 103) in 
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 20 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname = 'zangwil3';

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 100; -1; 2.5 ];             % x0
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'QUR2-AN-3-0';                % class

case 'cpsstr'

   cpsstr.name  = pname;
   cpsstr.eldom = { [ 1 2 3 ] };
   cpsstr.param = {};
   varargout{1} = cpsstr;

case 'objf'   % varargout = [ f, g, H ]

   x  = varargin{1};
   varargout = opm_eval_cpsf( pname, 'elobjf', x, { [ 1 2 3 ] }, nargout, {} );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   x  = varargin{2};
   r1 =  x(1)-x(2)+x(3);
   r2 = -x(1)+x(2)+x(3);
   r3 =  x(1)+x(2)-x(3);
   varargout{1} = r1^2 + r2^2 + r3^2;
   if ( nargout > 1 )
      J1 = [  1; -1;  1 ];
      J2 = [ -1;  1;  1 ];
      J3 = [  1;  1; -1 ];
      varargout{2} = 2 * ( J1*r1 + J2*r2 + J3*r3 );
      if ( nargout > 2 )
	 varargout{3} = 2 * ( J1*J1.' + J2*J2.' + J3*J3.' );
      end
   end
end

return

end