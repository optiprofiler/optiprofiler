%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                S2MPJ library for Python
%
%   Performs the runtime actions specific to S2MPJ, irrespective of the problem at hand.
%
%   Programming: Ph. Toint (this version 12 VI 2024)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = s2mpjlib( action, varargin )

switch ( action )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Computes the effective index of name in List, or add name to List if not in there already.
%   Return the index of varargin{1} in the varargin{2} dictionary and varargout{3} = 1 if varargin{2}
%   has been enlarged or 0 if varargin{3} was already present in varargin{2} at the call.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'ii'

   %  [ idx, List, new ] <== ( name, List )

   name = varargin{1};
   List = varargin{2};
   if ( isKey( List, name ) )
      varargout{1} = List( name );
      if ( nargout > 2 )
         varargout{3} = 0;
      end
   else
      varargout{1} = numEntries( List ) + 1;
      List( name ) = varargout{1};
      if ( nargout > 2 )
         varargout{3} = 1;
      end
   end
   varargout{2} = List;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Get the index of a nonlinear variable.  This implies adding it to the variables' dictionary ix_
%   if it is a new one with bounds, start point and types defined by their default settings (it is
%   too late to define problem-specific values).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'nlx'

%  [ iv, ix_, pb ] <== ( vname, ix_, pb, getxnames, xlowdef, xuppdef, x0def )

   [ iv, ix_, newvar] = s2mpjlib( 'ii', varargin{1} ,varargin{2} );
   pb   = varargin{3};
   if( newvar )
      pb.n = pb.n + 1;
      if ( varargin{4} )
         pb.xnames{iv} = varargin{1};
      end
      if ( ~isempty( varargin{4} ) )
         pb.xlower(iv,1) = varargin{4};
      else
         pb.xlower(iv,1) = 0.0;
      end
      if ( ~isempty( varargin{5} ) )
         pb.xupper(iv,1) = varargin{5};
      else
         pb.xupper(iv,1) = +Inf;
      end
      if ( isfield( pb, 'xtype' ) )
         pb.xtype(iv) = 'r';
      end
      if ( ~isempty( varargin{6} ) )
         pb.x0(iv,1) = varargin{6};
      else
         pb.x0(iv,1) = 0;
      end
   end
   varargout{1} = iv;
   varargout{2} = ix_;
   varargout{3} = pb;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Convert the real fields of a struct given by varargin{1} to multi-precision with a number of digits
%   given by varargin{2} (using vpa from the Symbolic Math package).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'convert'

%   converted_struct <== ( struct, precision )

    %  Set the precision for variable-precision arithmetic.
    
    digits( varargin{2} );

    %  Initialize the output structure.
    
    varargout{1} = varargin{1};

    % Get the field names of the structure to convert.
    
    fields = fieldnames( varargin{1} );

    for i = 1:length( fields )
        field = fields{i};
        value = varargin{1}.(field);

        if ( isnumeric( value ) )
        
            %  Check if the value has non-integer elements and convert non-integer elements to variable-precision
            %  arithmetic.

            if ( any( mod( value, 1 ) ~= 0 ) )
                varargout{1}.(field) = vpa( value );
%            else
%                varargout{1}.(field) = value;
            end

        % Recursively convert nested structures
            
%        elseif ( isstruct( value ) )
%            structOut.(field) = s2mpjlib( 'convert', value, precision );
        end
    end


   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Perform the main high-level action requested from a S2MPJ problem file, using information computed
%   during the 'setup' call and passed to this function in the pbm struct.  The nonlinear functions
%   are evaluated at x = varargin{2}. See the documentation of S2MPJ for details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case { 'fx', 'fgx', 'fgHx', 'cx', 'cJx', 'cJHx', 'cIx', 'cIJx', 'cIJHx', 'fHxv', 'cJxv', 'cIJxv', ...
       'Lxy','Lgxy', 'LgHxy', 'LIxy','LIgxy', 'LIgHxy', 'LHxyv','LIHxyv'}

   %  varargout <==( action, pbm, varargin )

   pbm  = varargin{1};

   %  Evaluate globals parameters, if present (as indicated by the pbm.has_globs flags).

   if ( pbm.has_globs(1) )
      pbm = feval( pbm.name, 'e_globs', pbm );
   end
   if ( pbm.has_globs(2) )
      pbm = feval( pbm.name, 'g_globs', pbm );
   end

   %  Call the appropriate action.

   switch ( action )
   case { 'fx', 'fgx', 'fgHx' }                      % varargin = { pbm, x }
      if ( ( isfield( pbm,'objgrps' ) && ~isempty( pbm.objgrps ) ) || ...
           ( isfield( pbm, 'H' ) && ~isempty( pbm.H ) ) )
         [varargout{1:nargout}] = evalgrsum( 1, pbm.objgrps, varargin{2}, pbm );
      else
         disp( ' ERROR: no objective function!' )
      end
   case 'fHxv'                                          % varargin = { problem, x, v }
      if ( ( isfield( pbm,'objgrps' ) && ~isempty( pbm.objgrps ) ) || ...
           ( isfield( pbm, 'H' ) && ~isempty( pbm.H ) ) )
         varargout{1} = evalHJv( 'Hv', pbm.objgrps, varargin{2}, varargin{3}, [], pbm );
      else
         disp( ' ERROR: no objective function!' )
      end
   case { 'cx', 'cJx', 'cJHx' }                      % varargin = { pbm, x }
      if ( isfield( pbm, 'congrps' ) && ~isempty( pbm.congrps ) )  % Check there are constraints!
         [varargout{1:nargout}] = evalgrsum( 0, pbm.congrps, varargin{2}, pbm );
      else
         disp( ' ERROR: no constraint!' )
      end
   case { 'cIx', 'cIJx', 'cIJHx' }                % varargin = { problem, x, I }
      if ( isfield( pbm, 'congrps' ) && ~isempty( pbm.congrps ) )  % Check there are constraints!
         [varargout{1:nargout}] = evalgrsum( 0, pbm.congrps( varargin{3} ), varargin{2}, pbm );
      else
         disp( ' ERROR: no constraint!' )
      end
   case 'cJxv'                                          % varargin = { pbm, x, v }
      if ( isfield( pbm, 'congrps' ) && ~isempty( pbm.congrps ) )  % Check there are constraints!
         varargout{1} = evalHJv( 'Jv', pbm.congrps, varargin{2}, varargin{3}, [], pbm );
      else
         disp( ' ERROR: no constraint!' )
      end
   case 'cIJxv'                                        % varargin = { pbm, x, v, I }
      if ( isfield( pbm, 'congrps' ) && ~isempty( pbm.congrps ) )  % Check there are constraints!
         varargout{1} = evalHJv( 'Jv', pbm.congrps( varargin{4} ), varargin{2}, varargin{3}, [], pbm );
      else
         disp( ' ERROR: no constraint!' )
      end
   case { 'Lxy', 'Lgxy', 'LgHxy' }                % varargin = { pbm, x, y }
      if ( ~isfield( pbm, 'objgrps' ) )
         pbm.objgrps = [];
      end
      if ( ~isfield( pbm, 'congrps' ) )
         pbm.congrps = [];
      end
      [varargout{1:nargout}] = evalLx( action, pbm.objgrps, pbm.congrps, varargin{2}, varargin{3}, [], pbm );
   case 'LHxyv'                                        % varargin = { pbm, x, y, v }
      if ( ~isfield( pbm, 'objgrps' ) )
         pbm.objgrps = [];
      end
      if ( ~isfield( pbm, 'congrps' ) )
         pbm.congrps = [];
      end
      varargout{1} = evalLx( action, pbm.objgrps, pbm.congrps, varargin{2}, varargin{3}, varargin{4}, pbm );
   case { 'LIxy', 'LIgxy', 'LIgHxy' }          % varargin = { pbm, x, y, I }
      if ( ~isfield( pbm, 'objgrps' ) )
         pbm.objgrps = [];
      end
      if ( ~isfield( pbm, 'congrps' ) )
         pbm.congrps = [];
      end
      y            = varargin{3};
      I            = varargin{4};
      [varargout{1:nargout}] = evalLx( action, pbm.objgrps, pbm.congrps(I), varargin{2}, y(I), pbm );
   case 'LIHxyv'                                      % varargin = { pbmm, x, y, v, I }
      if ( ~isfield( pbm, 'objgrps' ) )
         pbm.objgrps = [];
      end
      if ( ~isfield( pbm, 'congrps' ) )
         pbm.congrps = [];
      end
      y            = varargin{3};
      I            = varargin{5};
      varargout{1} = evalLx( action, pbm.objgrps, pbm.congrps(I), varargin{2}, y(I), varargin{4}, pbm );
   end

end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = evalgrsum( isobj, glist, x, pbm )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Evaluate the value of a sum of groups (and, if requested, that of of its gradient and Hessian) at x,
%   given the problem data avaialble in the pbm struct.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

debug = 0;%D

%  See if reduced precision is requested.

if ( isfield( pbm, 'ndigs' ) )
   ndigs   = pbm.ndigs;
   redprec = 1;
else
   redprec = 0;
end

%  Initializations

n  = length( x );
if ( isobj )
   fx = 0;
else
   m  = length( glist );
   cx = zeros( m, 1 );
   ic = 0;
end
if ( nargout > 1 )
   if ( redprec )
      if ( isobj )
         gx = 0*x;
      else
         Jx = sparse( m, n )*(0*x(1));
      end
   else
      if ( isobj )
         gx = zeros( n,1 );
      else
         Jx = sparse( m, n );
      end
   end
   if ( nargout > 2 )
      if ( isobj )
         Hx = sparse( n, n );
      else
         Hx = {};
      end
   end
end

%  Check for the presence and size of a linear term

if ( isfield( pbm, 'A' ) && ~isempty( pbm.A ) )
   [ sA1, sA2 ] = size( pbm.A );
   has_A = 1;
else
   has_A = 0;
end

%  Evaluate the quadratic term, if any.

if ( isobj && isfield( pbm, 'H' ) && ~isempty( pbm.H ) )
   Htimesx = pbm.H * x;
   switch( nargout )
   case 1
      fx = fx + 0.5 * x.' * Htimesx;
   case 2
      gx = gx + Htimesx;
      fx = fx + 0.5 * x.' * Htimesx;
   case 3
      Htimesx = pbm.H * x;
      gx = gx + Htimesx;
      fx = fx + 0.5 * x.' * Htimesx;
      Hx = Hx + pbm.H;
   end
end

if( debug )
   if ( isobj )
      fprintf( ' fx(quadratic) = %g\n',  fx );%D
      gxp = gx'%D
   else
      fprintf( ' cx(quadratic) = %g\n',  cx );%D
      full(Jx)%D
   end
end

%  Loop on the groups list

for iig = 1:length( glist )

   ig = glist( iig );

   %  Find the group's scaling.
   
   if ( isfield( pbm, 'gscale' ) )
      if ( ig <= length(pbm.gscale) &&  abs( pbm.gscale(ig) ) > 1.0e-15 )
         gsc = pbm.gscale(ig);
      else
         gsc = 1;
      end
   else
      gsc = 1;
   end

   %  Evaluate the linear term, if any.

   if ( isfield( pbm, 'gconst' ) )
      fin = -pbm.gconst( ig );
   else
      fin = 0;
   end
   switch( nargout )
   case 1
      if ( has_A && ig <= sA1 )
         fin = fin + pbm.A( ig, 1:sA2 ) * x(1:sA2);
      end
   case { 2, 3 }
      if ( redprec )
%         gin = 0*x;
         gin = sparse( n, 1 ) *(0*x(1));
      else
%         gin = zeros(n,1);
         gin = sparse( n, 1 );
      end
      if ( has_A && ig <= sA1 )
         gin(1:sA2) = pbm.A( ig, 1:sA2 ).';
         fin        = fin + gin.' * x;
      end
   end
   if ( nargout > 2 )
      if ( redprec )
         Hin = sparse( n, n )*(0*x(1));
      else
         Hin = sparse( n, n );
      end
   end

if( debug )
   fprintf( ' ig = %d  fin(linear) = %g\n', ig, fin );%D
   ginp = gin'%D
end

   %  Loop on the group's elements.
   %
   %  The explicit sequential scalar assignment of gradient and Hessian components may be
   %  necessary in this section because some elements may have repeated elemental variables,
   %  and only one occurence of such variables would be assigned by the vector assignments
   %  if a vector assignment were used. The case where explicit weights are present is treated
   %  separately in order to improve speed when weights are not present.

   if ( isfield( pbm, 'grelt' ) && ig <= length( pbm.grelt ) && ~isempty( pbm.grelt{ ig } ) )
      for iiel = 1:length( pbm.grelt{ ig } )       %  loop on elements
         iel    = pbm.grelt{ ig }( iiel );         %  the element's index
         irange = pbm.elvar{ iel };                %  the element's elemental variables
         efname = pbm.elftype{ iel };              %  the element's ftype

         if ( isfield( pbm, 'grelw' ) && ig <= length( pbm.grelw ) && ~isempty( pbm.grelw{ ig } ) )
            has_weights = 1;
            wiel        = pbm.grelw{ ig }( iiel );
         else
            has_weights = 0;
         end
         switch ( nargout )

         % Only the value is requested.
         
         case 1
            fiel = feval( pbm.name, efname, x( irange ), iel );
            if ( has_weights )
               fin = fin + wiel * fiel;
            else
               fin = fin + fiel;
            end

         %  The value and its gradient are requested.
         
         case 2
            [ fiel, giel ] = feval( pbm.name, efname, x( irange ), iel );
            if ( has_weights )
               fin = fin +  wiel * fiel;
               for ir = 1:length( irange )
                  ii  = irange( ir );
                  gin( ii ) = gin( ii ) + wiel * giel( ir );
               end
            else
               fin    = fin + fiel;
               for ir = 1:length( irange )
                  ii  = irange( ir );
                  gin( ii ) = gin( ii ) + giel( ir );
               end
            end

         %  The value, its gradient and its Hessian are requested.
         
         case 3
            [ fiel, giel, Hiel ] = feval( pbm.name, efname, x( irange ), iel );
            if ( has_weights )
               fin = fin +  wiel * fiel;
               for ir = 1:length( irange )
                  ii  = irange( ir );
                  gin( ii ) = gin( ii ) + wiel * giel( ir );
                  for jr = 1:length( irange )
                     jj  = irange( jr );
                     Hin( ii, jj ) = Hin( ii, jj ) + wiel * Hiel( ir, jr );
                  end
               end
            else
               fin = fin + fiel;
               for ir = 1:length( irange )
                  ii  = irange( ir );
                  gin( ii ) = gin( ii ) + giel( ir );
                  for jr = 1:length( irange )
                     jj  = irange( jr );
                     Hin( ii, jj ) = Hin( ii, jj ) + Hiel( ir, jr );
                  end
               end
            end
         end
      end
   end

if( debug )
   if ( isobj )
      fprintf( ' ig = %d  fin(nonlinear) = %g    fx = %g \n', ig, fin, fx );
      ginp = gin'
      gxp = gx'
   else
      fprintf( ' ig = %d  fin(nonlinear) = %g    cx = %g \n', ig, fin, cx );
      ginp = gin'
      full(Jx)
   end
end

   %  Evaluate the group function.
   %  1) the non-TRIVIAL case
   
   if ( isfield( pbm, 'grftype' )    && ig <= length( pbm.grftype ) && ...
       ~isempty( pbm.grftype{ ig } ) && ~strcmp( pbm.grftype{ ig}, 'TRIVIAL') )
      egname = pbm.grftype{ ig };            %  the group's ftype
      if ( isobj )
         switch ( nargout )
         case 1
            fx = fx + feval( pbm.name, egname, fin, ig ) / gsc;
         case 2
            [ fa, grada ] = feval( pbm.name, egname, fin, ig );
            fx = fx + fa / gsc;
            gx = gx + grada * gin /gsc ;
         case 3
           [ fa, grada, Hessa ] = feval( pbm.name, egname, fin, ig );
            fx   = fx + fa / gsc;
            gx   = gx + grada * gin / gsc;
            Hx   = Hx + ( ( Hessa * gin  ) * gin.' + grada * Hin ) / gsc;
         end
      else
         ic = ic + 1;
         switch ( nargout )
         case 1
            cx( ic ) = feval( pbm.name, egname, fin, ig ) / gsc;
         case 2
            [ fa, grada ] = feval( pbm.name, egname, fin, ig );
            cx( ic )      = fa / gsc;
            Jx(ic,1:n)    = grada * gin.' / gsc;
         case 3
            [ fa, grada, Hessa ] = feval( pbm.name, egname, fin, ig );
            cx( ic )   = fa / gsc;
            Jx(ic,1:n) = grada * gin.' / gsc;
            Hx{end+1}  = ( ( Hessa * gin ) * gin.' + grada * Hin ) / gsc;
         end
      end

   %  2) the TRIVIAL case: the group function is the identity
   
   else
      if ( isobj )
         switch ( nargout )
         case 1
            fx = fx + fin / gsc;
         case 2
            fx = fx + fin / gsc;
            gx = gx + gin / gsc;
         case 3
            fx = fx + fin / gsc;
            gx = gx + gin / gsc;
            Hx = Hx + Hin / gsc;
         end
      else
         ic = ic + 1;
         switch ( nargout )
         case 1
            cx( ic )   = fin  / gsc;
         case 2
            cx( ic )   = fin  / gsc;
            Jx(ic,1:n) = gin.'/ gsc;
         case 3
            cx( ic )   = fin  / gsc;
            Jx(ic,1:n) = gin.'/ gsc;
            Hx{end+1}  = Hin  / gsc;
         end
      end
   end

   if( debug )
      if( isobj)
         fprintf( ' ig = %d  fx(final) = %g\n', ig, fx);
%         gxp(final) = gx'%D
      else
         fprintf( ' ig = %d  cx(final) = %g\n', ig, cx);
         full(Jx)%D
      end
   end

end

%   Massage output.

if ( isobj )
   varargout{ 1 } = fx;
else
   varargout{ 1 } = cx;
end
if ( nargout > 1 )
   if ( isobj )
      varargout{ 2 } = gx;
   else
      varargout{ 2 } = Jx;
   end
   if ( nargout > 2 )
      varargout{ 3 } = Hx;
   end	 
end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function HJv = evalHJv( mode, glist, x, v, y, pbm )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Depending on mode:
%   mode = "Hv"  : evaluate the product of the objective's Hessian times v (glist = obj groups)
%   mode = "HIv" : evaluate the product of the constraints' Hessian  times v time the multiplier y
%                 (glist = cons groups)
%   mode = "Jv"  : evaluate the product of the constraints' Jacobian times v (glist = cons groups)
%   The vector y is unused (and unreferenced) for modes "Hv" and "Jv".
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

debug = 0;%D

%  See if reduced precision is requested.

if ( isfield( pbm, 'ndigs' ) )
   ndigs   = pbm.ndigs;
   redprec = 1;
else
   redprec = 0;
end

%  Initializations

n = length( x );
switch( mode )
case 'Hv'
   HJv = zeros( n, 1 );
case 'HIv'
   HJv = zeros( n, 1 );
   ic  = 0;
case 'Jv'
   HJv = zeros( length(glist), 1 );
   ic  = 0;
otherwise
   disp( [' ERROR: unknown mode', mode, ' at entry of evalHJv.' ] )
end

%  Check for the presence and size of a linear term.

if ( ~isempty( pbm.A ) )
   [ sA1, sA2 ] = size( pbm.A );
   has_A = 1;
else
   has_A = 0;
end

%  Evaluate the quadratic term, if any.

if ( strcmp( mode, 'Hv' ) && isfield( pbm, 'H' ) && ~isempty( pbm.H ) )
   HJv = HJv + pbm.H * v;
end
   
if( debug )
   fprintf( 'HJv(quadratic) = \n');%D
   HJv
end

%  Loop on the groups list

for iig = 1:length( glist )

   ig = glist( iig );

   %  Find the group's scaling.
   
   if ( isfield( pbm, 'gscale' ) )
      if ( ig <= length( pbm.gscale ) &&  abs( pbm.gscale( ig ) ) > 1.0e-15 )
         gsc = pbm.gscale(ig);
      else
         gsc = 1;
      end
   else
      gsc = 1;
   end

   %  Evaluate the linear term, if any.

   if ( isfield( pbm, 'gconst' ) )
      fin = -pbm.gconst( ig );
   else
      fin = 0;
   end
   if ( redprec )
      gin  = sparse( n, 1 )*(0*x(1));
   else
      gin  = sparse( n, 1 );
   end
   if ( has_A && ig <= sA1 )
      fin        = fin + pbm.A( ig, 1:sA2 ) * x(1:sA2);
      gin(1:sA2) = pbm.A( ig, 1:sA2 ).';
   end

if(debug )
   fprintf( ' ig = %d  fin(linear) = %g\n', ig, fin );%D
end
   %  Initialize the Hessian.

   if ( redprec )
      Hinv = 0*x;
   else
      Hinv = zeros( n, 1 );
   end

   %  Loop on the group's elements.
   %
   %  The explicit sequential scalar assignment of gradient and Hessian components may be
   %  necessary in this section because some elements may have repeated elemental variables,
   %  and only one occurence of such variables would be assigned by the vector assignments.

   if ( isfield( pbm, 'grelt' ) && ig <= length( pbm.grelt ) && ~isempty( pbm.grelt{ ig } ) )
      for iiel = 1:length( pbm.grelt{ ig } )
         iel    = pbm.grelt{ ig }( iiel );       %  the element's index
         irange = pbm.elvar{ iel };              %  the element's list of elemental variables
         efname = pbm.elftype{ iel };            %  the element's ftype
         if ( isfield( pbm, 'grelw' ) && ig <= length( pbm.grelw ) && ~isempty( pbm.grelw{ ig } ) )
            has_weights = 1;
            wiel        = pbm.grelw{ ig }( iiel );
         else
            has_weights = 0;
         end

         %  The group is an objective group.
         
         if ( ismember( mode, { 'Hv', 'HIv' } ) )
            [ fiel, giel, Hiel ] = feval( pbm.name, efname, x( irange ), iel );
            if ( has_weights )
               fin = fin + wiel * fiel;
               for ir = 1:length( irange )
                  ii  = irange( ir );
                  gin( ii ) = gin( ii )  + wiel * giel( ir );
                  for jr = 1:length(irange)
                     jj  = irange( jr );
                     Hinv( ii ) = Hinv( ii ) + wiel * Hiel( ir, jr ) * v( jj );
                  end
                end
            else
               fin = fin + fiel;
               for ir = 1:length( irange )
                  ii  = irange( ir );
                  gin( ii ) = gin( ii )  + giel( ir );
                  for jr = 1:length( irange )
                     jj  = irange( jr );
                     Hinv( ii ) = Hinv( ii ) + Hiel( ir, jr ) * v( jj );
                  end
               end
            end

         %  The groups is a constraint's group.
         
         else
            [ fiel, giel ] = feval( pbm.name, efname, x( irange ), iel );
            if ( has_weights )
               fin = fin + wiel * fiel;
               for ir = 1:length( irange )
                  ii  = irange( ir );
                  gin( ii ) = gin( ii ) + wiel * giel( ir );
               end
            else
               fin = fin + fiel;
               for ir = 1:length( irange )
                  ii = irange( ir );
                  gin( ii ) = gin( ii ) + giel( ir );
               end
            end
         end
      end
   end

   if( debug )
      fprintf( ' ig = %d Hinv(nonlinear) = \n', ig)
      Hinv
      fprintf( 'HJv = \n');%D
      HJv
   end

   %  Include contribution from the group function.
   
   if ( isfield( pbm, 'grftype' ) && ig <= length( pbm.grftype ) && ~isempty( pbm.grftype{ ig } ) )
      egname = pbm.grftype{ ig };          % the group's ftype
   else
      egname = 'TRIVIAL';
   end

   if ( ismember( mode, {'Hv', 'HIv' } ) )
      if ( strcmp( egname, 'TRIVIAL' ) )
         if ( strcmp( mode, 'HIv' ) )
            ic = ic + 1;
            HJv = HJv + y(ic) * Hinv / gsc;
         else
            HJv = HJv + Hinv / gsc;
         end
      else
         [ ~, grada, Hessa ] = feval( pbm.name, egname, fin, ig );
         if ( strcmp( mode, 'HIv' ) )
            ic = ic + 1;
            HJv  = HJv + y(ic) * ( ( Hessa * gin ) * ( gin.' * v ) + grada * Hinv ) / gsc;
         else
            HJv  = HJv + ( ( Hessa * gin ) * ( gin.' * v ) + grada * Hinv ) / gsc;
         end
      end
   else
      ic = ic + 1;
      if ( strcmp( egname, 'TRIVIAL' ) )
         HJv(ic) = gin.' * v / gsc;
      else
         [ ~, grada ] = feval( pbm.name, egname, fin, ig );
         HJv(ic)      = grada * gin.' * v / gsc;
      end
   end
end

if( debug )
   fprintf( ' ig = %d  HJv(final) = %g\n', ig, HJv);%D
end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = evalLx( action, gobjlist, gconlist, x, y, v, pbm )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Evaluate the value and derivatives of the Lagrangian function, or the product of the Lagrangian's
%   Hessian times a vector v, depending on "action".
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  Evaluate the objective function's part.

n = length(x);

if ( isempty( gobjlist ) && ~isfield( pbm, 'H' ) )
   switch ( action )
   case { 'Lxy', 'Lgxy', 'LgHxy', 'LIxy', 'LIgxy', 'LIgHxy' }
      varargout{1} = 0;
      if ( nargout > 1 )
         varargout{2} = zeros( n, 1 );
         if ( nargout > 2 )
            varargout{3} = sparse( n, n );
         end
      end
   case { 'LHxyv', 'LIHxyv' }
      varargout{1} = zeros( n, 1 );
   end
else
   switch ( action )
   case { 'Lxy', 'Lgxy', 'LgHxy', 'LIxy', 'LIgxy', 'LIgHxy' }
      [varargout{1:nargout}] = evalgrsum( 1, gobjlist, x, pbm );
   case { 'LHxyv', 'LIHxyv' }
      [varargout{1:nargout}] = evalHJv( 'Hv', gobjlist, x, v, [], pbm );
   end
end

%  Loop on the constraint's list.

switch( nargout )
case 1
   switch( action )
   case { 'Lxy', 'Lgxy', 'LgHxy', 'LIxy', 'LIgxy', 'LIgHxy' }
      varargout{1}    = varargout{1} + y' * evalgrsum( 0, gconlist, x, pbm );
   case { 'LHxyv', 'LIHxyv' }
      varargout{1} = varargout{1} + evalHJv( 'HIv', gconlist, x, v, y, pbm );
%      for ig = 1:length( gconlist )
%         varargout{1} = varargout{1} + y( ig ) * evalHJv( 'HIv', gconlist(ig), x, v, pbm );
%      end
   end
case 2
   [ coutvals{1:2} ]  = evalgrsum( 0, gconlist, x, pbm );
   varargout{ 1 }     = varargout{ 1 } + y' * coutvals{ 1 };
   varargout{ 2 }     = varargout{ 2 } + coutvals{ 2 }.' * y;
case 3
   [coutvals{1:3}]    = evalgrsum( 0, gconlist, x, pbm );
   varargout{ 1 }     = varargout{ 1 } + y' * coutvals{ 1 };
   varargout{ 2 }     = varargout{ 2 } +  coutvals{ 2 }.' * y;
   for  ig = 1:length( gconlist )
      varargout{ 3 }  = varargout{ 3 } + y( ig ) * coutvals{ 3 }{ ig };
   end
end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
