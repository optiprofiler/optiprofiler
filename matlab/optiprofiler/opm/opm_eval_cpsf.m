function fgH = opm_eval_cpsf( pname, fname, x, eldom, nout, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Evaluate the value of a CPS function (and, if requested, that of of its gradient and Hessian) at x,
%   given the function name and the list of its elemental variables for each of its elements.
%
%   INPUT:
%
%   pname   : the name of the considered OPM problem
%   fname   : the name of the considered CPS function within that problem
%   x       : the vector at which it must be evaluated
%   eldom   : a cell array whose i-th component is the list of elemental variables for element i
%   nout    : an output request:
%             if nout = 1, only the function value is requested
%             if nout = 2, the function value and its gradient are requested
%             if nout = 3, the function value, its gradient and its Hessian are requested
%   varargin: possible parameters.
%   OUTPUT:
%
%   fgH     : a cell array of length nout, with
%             fgH{1}: the function value at x
%             fgH{2}: the gradient value at x (if nout > 1)
%             fgH{3}: the Hessian  value at x (if nout > 2), as a sparse matrix.
%
%   Programming : Ph. Toint and S. Gratton, July 2018
%   This version: 21 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n       = length( x );
f       = 0;
if ( nout > 1 )
   g    = zeros( n, 1 );
   if ( nout > 2 )
      H = sparse( n, n );
   end
end

%  Loop on the function's elements.

for iel = 1:length( eldom )
   irange = eldom{ iel };
   switch ( nout )
   case 1
      fiel = feval( pname, fname, iel, x( irange ), varargin{:} );
      f  = f + fiel;
   case 2
      [ fiel, giel ] = feval( pname, fname, iel, x( irange ), varargin{:} );
      f  = f + fiel;
      g( irange ) = g( irange ) + giel;
  case 3
      [ fiel, giel, Hiel ] = feval( pname, fname, iel, x( irange ), varargin{:} );
      f  = f + fiel;
      g( irange ) = g( irange ) + giel;
      H( irange, irange ) = H( irange, irange ) + Hiel;
   end
end

%   Massage output.

fgH      = cell( nout, 1 );
fgH{ 1 } = f;
if ( nout > 1 )
   fgH{ 2 } = g;
   if ( nout > 2 )
      fgH{ 3 } = H;
   end	 
end

return
