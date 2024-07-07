%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = himm33( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   A nonlinear least-squares problem in three variables.
%
%   Source: problem 87 (p. 67) in
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 21 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'himm33';

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 0.5; 0.5 ];                 % x0
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

   x     = varargin{1};
   varargout = opm_eval_cpsf( pname, 'elobjf', x, { [ 1 2 ] }, nargout, {} );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   x    = varargin{2};
   ee   = exp(-(x(1)+x(2)));
   rr   = ( 2*x(1)^2 + 3*x(2)^2 );
   varargout{1} = ee*rr;
   if ( nargout > 1 )
      Jee = -ee * [ 1; 1 ];
      Jrr = [ 4*x(1); 6*x(2) ];
      varargout{2} = Jee*rr+ee*Jrr;
      if ( nargout > 2 )
         Hee = ee* [ 1, 1;
	             1, 1 ];
         Hrr = [ 4, 0;
	         0, 6 ];
         varargout{3} = Jee*Jrr.'+rr*Hee+Jrr*Jee.'+ee*Hrr;
      end
   end
end

return

end
