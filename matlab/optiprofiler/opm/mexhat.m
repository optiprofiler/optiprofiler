%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = mexhat( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The mexican hat problem with penalty parameter 0.00001
%
%   Source:
%   A.A. Brown and M. Bartholomew-Biggs,
%      "Some effective methods for unconstrained optimization based on
%      the solution of ordinary differential equations",
%      Technical Report 178, Numerical Optimization Centre, Hatfield
%      Polytechnic, (Hatfield, UK), 1987.
%
%   Ph. Toint 28 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname = 'mexhat';

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 0.86, 0.72 ]';              % x0
   varargout{2} = [-1.1171526 -0.0898793 ];     % fstar
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

   varargout = opm_eval_cpsf( pname, 'elobjf', varargin{1}, { [1:2] }, nargout, {} );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   x  = varargin{2};
   r1 = (x(1)-1)^2 + (x(2)-1)^2;
   r2 = sqrt(1e5)*(-0.02 + x(2)-x(1)^2);
   varargout{1} = r1^2 + r2^2;
   if( nargout > 1 )
     J1 = [2*(x(1)-1);2*(x(2)-1)];
     J2 = sqrt(1e5)*[ -2*x(1); 1 ];
     varargout{2} = 2 * ( J1*r1 + J2*r2 );
     if ( nargout > 2 )
        H1 = [ 2, 0; 0, 2 ];
	H2 = sqrt(1e5)*[ -2, 0; 0, 0 ];
	varargout{3}= 2*( J1*J1.' + r1*H1 + J2*J2.' + r2*H2 );
     end
   end
end

return

end