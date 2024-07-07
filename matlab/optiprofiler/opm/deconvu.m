%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = deconvu( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The DECONVU problem
%
%   A problem arising in deconvolution analysis.
%
%   Source:
%     J.P. Rasson, private communication, 1986.
%
%   Ph. Toint 21 VII 2021.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'deconvu';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   ssg = [ 0.01 0.02 0.4 0.6 0.8 3.0 0.8 0.6 0.44 0.01 0.01 ];
   varargout{1} = [ zeros(40,1); ssg' ];          % x0
   varargout{2} = NaN;                            % fstar
   varargout{3} = '';                             % xtype
   varargout{4} = [];                             % xlower
   varargout{5} = [];                             % xupper
   varargout{6} = [];                             % clower
   varargout{7} = [];                             % cupper
   varargout{8} = 'SUR2-MN-51-0';                  % class

case 'cpsstr'

   for iel = 1:40
      isg          = [ 40+1:40+min(11,iel) ];
      ic           = [ iel:-1:iel-length(isg)+1 ];
      eldom{ iel } = [ isg ic ];
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

   tr = [ 0.0      0.0    1.6e-03  5.4e-03 7.02e-02     0.1876 0.332 0.764 0.932  0.812  ...
          0.3464   0.2064 8.3e-02  3.4e-02 6.179999e-02 1.2    1.8   2.4   9.0     2.4   ...
	  1.801    1.325  7.62e-02 0.2104  0.268        0.552  0.996 0.36  0 .24   0.151 ...
	  2.48e-02 0.2432 0.3602   0.48    1.8          0.48   0.36  0.264 6.0e-03 6.0e-03 ];

   iel = varargin{1};
   x   = varargin{2};
   nx  = length(x);
   nt  = nx / 2;
   riel   = -tr( iel );
   for i = 1:nt
      riel = riel + x(i)*x(nt+i);
   end
   varargout{1} = riel^2;
   if ( nargout > 1 )
      Jiel = zeros( nx, 1 );
      for i = 1:nt
         Jiel(i)    = Jiel(i)    + x( nt+i );
	 Jiel(nt+i) = Jiel(nt+i) + x( i );
      end
      varargout{2} = 2* Jiel * riel;
      if ( nargout > 2 )
         Hiel = sparse( nx, nx );
      for i = 1:nt
	 Hiel( i, nt+i ) = Hiel( i, nt+i ) + 1;
	 Hiel( nt+i, i ) = Hiel( nt+i, i ) + 1;
      end
      varargout{3} = 2 * ( Jiel*Jiel' + riel*Hiel );
   end
end

return

end