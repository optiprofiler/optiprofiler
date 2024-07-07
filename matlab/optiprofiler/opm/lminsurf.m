%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = lminsurf( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   THe linear minimum surface problem.
%   This problem comes from the discretization of the minimum surface
%   problem on the unit square: given a set of boundary conditions on
%   the four sides of the square, one must find the surface which
%   meets these boundary conditions and is of minimum area.
%
%   The unit square is discretized into (p-1)^2 little squares. The
%   heights of the considered surface above the corners of these little
%   squares are the problem variables,  There are p^2 of them.
%   Given these heights, the area above a little square is
%   approximated by the
%     S(i,j) = sqrt( 1 + 0.5(p-1)^2 ( a(i,j)**2 + b(i,j)^2 ) ) / (p-1)^2
%   where
%     a(i,j) = x(i,j) - x(i+1,j+1)
%   and
%     b(i,j) = x(i+1,j) - x(i,j+1)
%
%   In the Linear Mininum Surface, the boundary conditions are given
%   as the heights of a given plane above the square boundaries.  This
%   plane is specified by its height above the (x,y,0) plane, and its slopes
%   along the first and second coordinate directions in the plane.
%
%   The number of variables must be a square and be at least 9.
%
%   Source:
%      A Griewank and Ph. Toint,
%      "Partitioned variable metric updates for large structured 
%      optimization problems",
%      Numerische Mathematik 39:429-448, 1982.
%
%   If the dimension is unspecified, the default n = 16 is chosen.
%
%   Ph. Toint 22 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'lminsurf';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class ]

   if ( length( varargin ) )
      n     = varargin{1};
      nsq   = sqrt( n );
      error = ( n < 9 || nsq ~= round( nsq ) );
      if ( error )
         disp( [ ' ERROR in lminsurf: n = ', int2str(n), ' but should be a square and at least 9!' ] )
      end
   else
      n   = 16;
      nsq = 4;
   end
   nsqm1      = nsq - 1;
   h          = 1 / nsqm1;
   x0         = zeros( n, 1 );
   xlower     = -Inf * ones( n, 1 );
   xupper     =  Inf * ones( n, 1 );
   for iy = 1:nsq
      if ( iy == 1 )
         for ix = 1:nsq
            t            = (ix-1)*h;
            x0( ix )     = 1 + 8*t;
	    xlower( ix ) = x0( ix );
	    xupper( ix ) = x0( ix );
         end
      elseif ( iy == nsq )
         for ix = 1:nsq
            ii           = ix + (nsq-1)*nsq;
	    t            = (ix-1)*h;
            x0( ii )     = 5 + 8*t;
	    xlower( ii ) = x0( ii );
	    xupper( ii ) = x0( ii );
         end
      else
	 ii           = (iy-1)*nsq+1;
	 t            = (iy-1)*h;
         x0( ii )     = 1 + 4*t;
	 xlower( ii ) = x0( ii );
	 xupper( ii ) = x0( ii );
	 ii           = iy*nsq;
         x0( ii )     = 9 + 4*t;
	 xlower( ii ) = x0( ii );
	 xupper( ii ) = x0( ii );
      end
   end
   varargout{1} = x0;
   varargout{2} = 9;              % fstar
   varargout{3} = '';             % xtype
   varargout{4} = xlower;
   varargout{5} = xupper;
   varargout{6} = [];             % clower
   varargout{7} = [];             % cupper
   varargout{8} = 'OXR2-MY-V-0';  % class

case 'cpsstr'

   nsq   = sqrt( varargin{1} );
   nsqm1 = nsq - 1;
   eldom = cell( nsqm1^2, 1 );
   for ii = 1:nsqm1
      for jj = 1:nsqm1
         iel          = (ii-1)*nsqm1 + jj;
         ixl          = (ii-1)*nsq   + jj;
	 eldom{ iel } = [ ixl ixl+1 ixl+nsq ixl+nsq+1 ];
      end
   end
   cpsstr.name  = pname;
   cpsstr.eldom = eldom;
   varargout{1} = cpsstr;

case 'objf'   % varargout = [ f, g, H ]

   x = varargin{1};
   n = length( x );
   if ( nargin > 2 && isstruct( varargin{2} ) && isfield( varargin{2}, 'eldom' ) );
      cpsstr = varargin{2};
   else
      cpsstr = problem( 'cpsstr', n );
   end
   varargout = opm_eval_cpsf( pname, 'elobjf', x, cpsstr.eldom, nargout, n );

   %   Zero derivatives of fixed variables.

   nsq = sqrt( n );
   if ( nargout > 1 )
      for iy = 1:nsq
         if ( iy == 1 )
            for ix = 1:nsq
	       varargout{2}(ix) = 0;
               if ( nargout > 2 )
	          varargout{3}( ix, 1:n ) = 0;
	          varargout{3}( 1:n, ix ) = 0;
	       end
            end
         elseif ( iy == nsq )
            for ix = 1:nsq
               ii           = ix + (nsq-1)*nsq;
	       varargout{2}(ii) = 0;
               if ( nargout > 2 )
	          varargout{3}( ii, 1:n ) = 0;
	          varargout{3}( 1:n, ii ) = 0;
	       end
            end
         else
	    ii           = (iy-1)*nsq+1;
	    varargout{2}(ii) = 0;
            if ( nargout > 2 )
	       varargout{3}( ii, 1:n ) = 0;
	       varargout{3}( 1:n, ii ) = 0;
	    end
  	    ii           = iy*nsq;
	    varargout{2}(ii) = 0;
            if ( nargout > 2 )
	       varargout{3}( ii, 1:n ) = 0;
	       varargout{3}( 1:n, ii ) = 0;
	    end
         end
      end
   end

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   x     = varargin{2};
   n     = varargin{3};
   nsq   = sqrt( n );
   nsqm1 = nsq - 1;
   nel   = nsqm1^2;
   ri    = 1 + 0.5 * nel * ( ( x(1) - x(4) )^2 + ( x(2) - x(3) )^2 );
   fiel  =  sqrt( ri );
   varargout{1}  = fiel/nel;
   if ( nargout > 1 )
      gri  = [ x(1)-x(4); x(2)-x(3); x(3)-x(2); x(4)-x(1) ];
      varargout{2} = ( 0.5 / fiel ) * gri;
      if ( nargout > 2 )
         varargout{3} = ( 0.5 / fiel ) * [  1,  0,  0, -1;
                                                  0,  1, -1,  0;
	  		                          0, -1,  1,  0;
	                                         -1,  0,  0,  1 ]  + ...
                        (-0.25 * nel / (ri*fiel) ) * gri * gri.';
%         varargout{3} = varargout{3}/nel;
      end
   end
   
end

return

end