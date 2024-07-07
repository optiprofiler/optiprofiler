%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = trigger( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The trigger circuit problem by Rheinboldt, as a function
%   of the input voltage.  This is a nonlinear least-squares.
%
%   Source:
%      G. Poenisch and H. Schwetlick,
%      "Computing Turning Points of Curves Implicitly Defined by
%      Nonlinear Equations Depending on a Parameter",
%      Computing 26:107-121, 1981.
%   Also problem 12 in
%   J.J. More',
%      "A collection of nonlinear model problems"
%      Proceedings of the AMS-SIAM Summer seminar on the Computational
%      Solution of Nonlinear Systems of Equations, Colorado, 1988.
%      Argonne National Laboratory MCS-P60-0289, 1989.
%   but there are several typos in that source.
%
%
%   Ph. Toint 27 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'trigger';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   invalue      = 0.322866124;
   varargout{1} = [ invalue; 0.2; 0.6; 0.2; 0.2; 0.6; 9.6];   % x0
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [ invalue; -Inf*ones(6,1) ];  % xlower
   varargout{5} = [ invalue;  Inf*ones(6,1) ];  % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SXR2-AN-7-0';                % class

case 'cpsstr'

   eldom{1}     = [ 1 2 3 7 ];
   eldom{2}     = [ 1 2 6 ];
   eldom{3}     = [ 1 3 4 ];
   eldom{4}     = [ 3 4 5 ];
   eldom{5}     = [ 4 5 6 ];
   eldom{6}     = [ 2 5 6 1 3 ];
   R  = [ 10000. 39. 51. 10. 25.5 1. 0.62 13. 0.201 ];
   A  = zeros( 6, 6 );
   A(1,1) = 1/R(1)+1/R(2)+1/R(3);
   A(1,2) = 1/R(2) - 1;
   A(2,2) = 1/R(2);
   A(2,6) = 1/R(4)-1;
   A(3,1) = 1/R(1)-1;
   A(3,3) = 1/R(1)+1/R(5);
   A(3,4) = 1/R(5)-1;
   A(4,4) = 1/R(5)+1/R(6)+1/R(7);
   A(4,5) = 1/R(6)-1;
   A(5,5) = 1/R(6)+1/R(8);
   A(5,6) = 1/R(8)-1;
   A(6,6) = 1/R(4)+1/R(8)+1/R(9);
   cpsstr.name  = pname;
   cpsstr.eldom = eldom;
   cpsstr.param = { A, R };
   varargout{1} = cpsstr;

case 'objf'   % varargout = [ f, g, H ]

   x = varargin{1};
   if ( nargin > 2 && isstruct( varargin{2} ) && isfield( varargin{2}, 'eldom' ) )
      cpsstr = varargin{2};
   else
      cpsstr = problem( 'cpsstr', length( x ) );
   end
   varargout = opm_eval_cpsf( pname, 'elobjf', x, cpsstr.eldom, nargout, cpsstr.param{:} );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   iel  = varargin{1};
   x    = varargin{2};
   A    = varargin{3};
   R    = varargin{4};
   b1   = 5.6e-8;
   b2   = 1962.0;
   switch ( iel )
   case 1
      riel         = A(1,1)*x(1) + A(1,2)*x(2) + A(3,1)*x(3) + x(4)/R(2);
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = [ A(1,1); A(1,2); A(3,1); 1/R(2) ];
	 varargout{2} = 2 * Jiel * riel;
         if ( nargout > 2 )
	    varargout{3} = 2* Jiel*Jiel.';
	 end
      end
   case 2
      e2           = b1*exp(25*(x(2)-1));
      riel         = A(1,2)*x(1) + A(2,2)*x(2) + A(2,6)*x(3) + e2;
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = [ A(1,2); A(2,2)+25*e2; A(2,6) ];
	 varargout{2} = 2 * Jiel * riel;
         if ( nargout > 2 )
	    Hiel = [ 0,   0,    0;
	             0, 625*e2, 0;
		     0,   0,    0 ];
	    varargout{3} = 2 * ( Jiel*Jiel.' +  riel*Hiel );
	 end
      end
   case 3
      riel = A(3,1)*x(1) + A(3,3)*x(2) + A(3,4)*x(3);
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = [ A(3,1); A(3,3); A(3,4) ];
	 varargout{2} = 2 * Jiel * riel;
         if ( nargout > 2 )
	    varargout{3} = 2* Jiel*Jiel.';
	 end
      end
   case 4
      riel         = A(3,4)*x(1) + A(4,4)*x(2) + A(4,5)*x(3);
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = [ A(3,4); A(4,4); A(4,5) ];
	 varargout{2} = 2 * Jiel * riel;
         if ( nargout > 2 )
	    varargout{3} = 2* Jiel*Jiel.';
	 end
      end
   case 5
      e2           = b1*exp(25*(x(2)-1));
      riel         = A(4,5)*x(1) + A(5,5)*x(2) + A(5,6)*x(3) + e2;
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = [ A(4,5); A(5,5)+25*e2; A(5,6) ];
	 varargout{2} = 2 * Jiel * riel;
         if ( nargout > 2 )
	    Hiel = [ 0,    0  , 0 ;
	             0, 625*e2, 0 ;
	 	     0,    0  , 0 ];
	    varargout{3} = 2 * ( Jiel*Jiel.' +  riel*Hiel );
	 end
      end
   case 6
      u    = x(5)-x(4);
      b2sq = b2^2;
      bu2  = 1+b2sq*u^2;
      riel = A(2,6)*x(1) + A(5,6)*x(2) + A(6,6)*x(3) + 7.65*atan(b2*u);
      varargout{1} = riel^2;
      if ( nargout > 1 )
         Jiel = [ A(2,6); A(5,6); A(6,6); 7.65*b2/bu2*[-1;1] ];
	 varargout{2} = 2 * Jiel * riel;
         if ( nargout > 2 )
	    Hiel = zeros( 5, 5 );
	    Hiel( [4 5],[4 5] ) = -15.3*b2sq*b2*u/bu2^2*[  1, -1;
		                                          -1,  1 ];
	    varargout{3} = 2 * ( Jiel*Jiel.' +  riel*Hiel );
	 end
      end
   end

end

return

end