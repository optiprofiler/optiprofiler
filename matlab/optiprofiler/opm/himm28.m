%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = himm28( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   A nonlinear least-squares problem in two variables.
%
%   Source: problem  6 (p. 63) in
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 21 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'himm28';

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 1; 1 ];                     % x0
   varargout{2} = 0;                            % fstar
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

   x     = varargin{1};
   varargout = opm_eval_cpsf( pname, 'elobjf', x, { [ 1 2 ] }, nargout, {} );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   x    = varargin{2};
   r1   = x(1)^2+x(2)-11;
   r2   = x(1)+x(2)^2-7;
   varargout{1} = r1^2 + r2^2;
   if ( nargout > 1 )
      J1 = [ 2*x(1); 1 ];
      J2 = [ 1; 2*x(2) ];
      varargout{2} = 2 * ( J1*r1 + J2*r2 );
      if ( nargout > 2 )
         H1 = [ 2, 0; 0, 0 ];
	 H2 = [ 0, 0; 0, 2 ];
         varargout{3} = 2*(J1*J1.'+r1*H1+J2*J2.'+r2*H2);
      end
   end
end

return

end
