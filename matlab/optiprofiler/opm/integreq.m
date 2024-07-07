%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = integreq( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The discrete integral problem.
%
%   Source: problem 29 in
%      J.J. More, B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%   Also problem 165 (p. 74) in
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   If the dimension is unspecified, the default n = 10 is chosen.
%
%   Ph. Toint 22 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'integreq';

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 2 )
         disp( [ ' ERROR in integreq: n = ', int2str(n), ' must be > 1!' ] )
      end
   else
      n = 10;
   end
   t            = [1:n].'/(n+1);
   varargout{1} = t.*(t-ones(n,1));             % x0
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AN-V-0';                % class

case 'cpsstr'

   cpsstr.name  = pname;
   cpsstr.eldom = { [ 1:varargin{1} ] };
   varargout{1} = cpsstr;
   
case 'objf'   % varargout = [ f, g, H ]

   x = varargin{1};
   n = length( x );
   varargout = opm_eval_cpsf( pname, 'elobjf', x, { [ 1:n ] }, nargout, n );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   x    = varargin{2};
   n    = varargin{3};
   varargout{1} = 0;
   if ( nargout > 1 )
      varargout{2} = zeros(n,1);
      if ( nargout > 2 )
         varargout{3} = zeros(n,n);
      end
   end
   t    = [1:n].'/(n+1);
   xtmp = x + t + ones(n,1);
   for i = 1:n
      sum1 = sum( t(1:i) .*xtmp(1:i).^3);
      sum2 = sum((ones(n-i,1)-t(i+1:n)).*xtmp(i+1:n).^3);
      ti   = t(i);
      ri   = x(i) + 0.5*((1-ti)*sum1 + ti*sum2);
      varargout{1} = varargout{1} + ri^2;
      if ( nargout > 1 )
	 Ji = [ 1.5*(1-ti)*t(1:i).*xtmp(1:i).^2;
	        1.5*ti*(ones(n-i,1)-t(i+1:n)).*xtmp(i+1:n).^2 ];
         Ji(i) = Ji(i)+1;
         varargout{2} = varargout{2} + 2 * Ji * ri;
         if ( nargout > 2 )
	    Hi = sparse( n, n );
	    Hi = Hi+3*diag([(1-ti)*t(1:i).*xtmp(1:i); ti*(ones(n-i,1)-t(i+1:n)).*xtmp(i+1:n)] );
	    varargout{3} = varargout{3} + 2 * ( Ji * Ji.' + ri * Hi );
	 end
      end
   end
end

return

end