%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = gulf( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The Gulf research and development function for m = 99. 
%
%   Source: problem 11 in
%      J.J. More', B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%   Also problem 27 (p. 57) in
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 20 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'gulf';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 5; 2.5; 0.15];               % x0
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AN-3-0';                % class

case 'cpsstr'

   for iel = 1:99
      eldom{ iel } = [ 1 2 3 ];
   end
   cpsstr.name  = pname;
   cpsstr.eldom = eldom;
   cpsstr.param = {};
   varargout{1} = cpsstr;

case 'objf'   % varargout = [ f, g, H ]

   x = varargin{1};
   if ( nargin > 2 && isstruct( varargin{2} ) && isfield( varargin{2}, 'eldom' ) );
      cpsstr = varargin{2};
   else
      cpsstr = problem( 'cpsstr', length( x ) );
   end
   varargout = opm_eval_cpsf( pname, 'elobjf', x, cpsstr.eldom, nargout, cpsstr.param{:} );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   iel   = varargin{1};
   x     = varargin{2};
   ymx2  = 25 + (-50*log(0.01*iel))^(2/3) - x(2);
   a     = abs( ymx2 )^(x(3)) / x(1);
   expma = exp( -a );
   riel  = expma - 0.01 * iel;
   varargout{1} = riel^2;
   if ( nargout > 1 )
      aexpma = a * expma;
      lnymx2 = log(abs(ymx2));
      Jiel = [ aexpma/x(1); x(3)*aexpma/ymx2; -aexpma*lnymx2];
      varargout{2} = 2 * Jiel * riel;
      if ( nargout > 2 )
         am1  = a - 1;
         aln  = a*lnymx2;
         Hiel = [ (a-2)*aexpma/x(1)^2, x(3)*am1*aexpma/(x(1)*ymx2), -aln*aexpma/x(1);
	           x(3)*am1*aexpma/(x(1)*ymx2),x(3)*aexpma*(1+x(3)*am1)/ymx2^2,aexpma*(1+x(3)*aln)/ymx2;
	          -aln*aexpma/x(1),aexpma*(1+x(3)*aln)/ymx2,aln*lnymx2*expma*am1 ];
	 varargout{3} = 2 * ( Jiel*Jiel.' + riel*Hiel );
      end
   end
   
end

return

end