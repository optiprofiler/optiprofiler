%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = vibrbeam( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   A nonlinear least-squares problem arising from laser-Doppler
%   measurements of a vibrating beam.  The data correspond to a simulated
%   experiment where two laser-Doppler velocimeters take measurements
%   at random points along the centreline of the beam.  These measurements
%   consist of a position (x), an incident angle (p) and the magnitude
%   of the velocity along the line of sight (v).
%   The problem is then to fit
%
%                         2      3                    2     3
%       v = (c + c x + c x  + c x ) cos[ d + d x + d x + d x  - p ]
%             0   1     2      3          0   1     2     3
%           <---- magnitude ----->       <------ phase ----->
%
%   in the least-squares sense.
%
%   Source: 
%   a modification of an exercize for L. Watson course on LANCELOT in
%   the Spring 1993. Compared to the original proposal, the unnecessary
%   elements were removed as well as an unnecessary constraint on the phase.
%
%   SIF input: Ph. L. Toint, May 1993, based on a proposal by
%              D. E. Montgomery, Virginia Tech., April 1993.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pname   = 'vibrbeam';
problem = str2func( pname );

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

%                    c0    c1  c2 c3 d0   d1 d2 d3
   varargout{1} = [ -3.5; 1.0; 0; 0; 1.7; 0; 0; 0 ]; % x0
   varargout{2} = '???';                          % fstar
   varargout{3} = '';                             % xtype
   varargout{4} = [];                             % xlower
   varargout{5} = [];                             % xupper
   varargout{6} = [];                             % clower
   varargout{7} = [];                             % cupper
   varargout{8} = 'SUR2-MN-8-0';                  % class

case 'cpsstr'

   nel = 30;
   for iel = 1:nel
      eldom{iel} = [1:8];
   end
%  Positions

   pos  = [ 39.1722;  53.9707; 47.9829; 12.5925; 16.5414;
            18.9548;  27.7168; 31.9201; 45.6830; 22.2524;
	    33.9805;   6.8425; 35.1677; 33.5682; 43.3659;
	    13.3835;  25.7273; 21.0230; 10.9755;  1.5323;
	    45.4416;  14.5431; 22.4313; 29.0144; 25.2675;
	    15.5095;   9.6297;  8.3009; 30.8694; 43.3299 ];

%  Velocity magnitude

   vel = [  -1.2026;   1.7053;  0.5410;  1.1477;  1.2447;
             0.9428;  -0.1360; -0.7542; -0.3396;  0.7057;
	    -0.8509;  -0.1201; -1.2193; -1.0448; -0.7723;
	     0.4342;   0.1154;  0.2868;  0.3558; -0.5090;
	    -0.0842;   0.6021;  0.1197; -0.1827;  0.1806;
	     0.5395;   0.2072;  0.1466; -0.2672; -0.3038 ];
	  
%  Angle of incidence

   ang = [   2.5736;   2.7078;  2.6613;  2.0374;  2.1553;
             2.2195;   2.4077;  2.4772;  2.6409;  2.2981;
	     2.5073;   1.8380;  2.5236;  2.5015;  2.6186;
	     0.4947;   0.6062;  0.5588;  0.4772;  0.4184;
	     0.9051;   0.5035;  0.5723;  0.6437;  0.6013;
	     0.5111;   0.4679;  0.4590;  0.6666;  0.8630 ];

   cpsstr.name  = pname;
   cpsstr.eldom = eldom;
   cpsstr.param = { pos, vel, ang };
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

   iel = varargin{1};
   x   = varargin{2};
   pos = varargin{3};
   vel = varargin{4};
   ang = varargin{5};
   riel = -vel(iel);
   switch ( nargout )
   case 1
      for j = 1:4
         fe   = vibrbeam( 'fun', [x(5:8);x(j)], pos(iel), ang(iel) );
         riel = riel + fe;
      end
      varargout{1} = riel^2;
   case 2
      Jiel = zeros(8,1);
      for j = 1:4
         [ fe, ge ] = vibrbeam( 'fun', [x(5:8);x(j)], pos(iel), ang(iel) );
         riel       = riel + fe;
         Jiel(5:8)  = Jiel(5:8) + ge(1:4);
	 Jiel(j)    = Jiel(j)   + ge(5);
      end
      varargout{1} = riel^2;
      varargout{2} = 2*Jiel*riel;
   case 3
      Jiel = zeros(8,1);
      Hiel = zeros(8,8);
      for j = 1:4
         [ fe, ge, He ] = vibrbeam( 'fun', [x(5:8);x(j)], pos(iel), ang(iel) );
         riel       = riel + fe;
         Jiel(5:8)  = Jiel(5:8) + ge(1:4);
	 Jiel(j)    = Jiel(j)   + ge(5);
	 Hiel(5:8,5:8) = Hiel(5:8,5:8) + He(1:4,1:4);
	 Hiel(5:8,j)   = Hiel(5:8,j) + He(1:4,5);
	 Hiel(j,5:8)   = Hiel(j,5:8) + He(5,1:4);
	 Hiel(j,j)     = Hiel(j,j)   + He(5,5);
      end
      varargout{1} = riel^2;
      varargout{2} = 2*Jiel*riel;
      varargout{3} = 2*( Jiel*Jiel.' + riel*Hiel);
   end

case 'fun'

   x = varargin{1};
   y = varargin{2};
   q = varargin{3};
   phi    =  x(1)+y*(x(2)+y*(x(3)+y*x(4)))-q;
   cosphi = cos( phi );
   bcos   = x(5) * cosphi;
   varargout{1} = bcos;
   if ( nargout > 1 )
      y2 = y * y;
      y3 = y * y2;
      y4 = y2 * y2;
      y5 = y2 * y3;
      y6 = y3 * y3;
      sinphi = sin( phi );
      bsin   = x(5) * sinphi;
      varargout{2} = [-bsin; - bsin*y; - bsin*y2; -bsin*y3; cosphi ];
      if ( nargout > 2 )
         varargout{3} = [ -bcos     -bcos*y   -bcos*y2   -bcos*y3   -sinphi;
	                  -bcos*y   -bcos*y2  -bcos*y3   -bcos*y4   -sinphi*y;
			  -bcos*y2  -bcos*y3  -bcos*y4   -bcos*y5   -sinphi*y2;
			  -bcos*y3  -bcos*y4  -bcos*y5   -bcos*y6   -sinphi*y3;
			  -sinphi   -sinphi*y -sinphi*y2 -sinphi*y3      0     ];
       end
    end
end

return

end
