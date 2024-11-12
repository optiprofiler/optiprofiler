%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                S2MPJ library for Matlab
%
%   Performs the runtime actions specific to S2MPJ, irrespective of the problem at hand.
%
%   Programming: Ph. Toint (this version 9 XI 2024)
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
      switch( class( List ) )
      case 'dictionary'
         varargout{1} = numEntries( List ) + 1;
      case 'containers.Map'
         varargout{1} = List.Count + 1;
      end      
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

case { 'fx', 'fgx', 'fgHx', 'cx', 'cJx', 'cJHx', 'cIx', 'cIJx', 'cIJHx', 'fHxv', 'cJxv', 'cJtxv', ...
       'cIJxv', 'cIJtxv', 'Lxy','Lgxy', 'LgHxy', 'LIxy','LIgxy', 'LIgHxy', 'LHxyv','LIHxyv'}

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
   case 'cJtxv'                                          % varargin = { pbm, x, v }
      if ( isfield( pbm, 'congrps' ) && ~isempty( pbm.congrps ) )  % Check there are constraints!
         varargout{1} = evalHJv( 'Jtv', pbm.congrps, varargin{2}, varargin{3}, [], pbm );
      else
         disp( ' ERROR: no constraint!' )
      end
   case 'cIJxv'                                        % varargin = { pbm, x, v, I }
      if ( isfield( pbm, 'congrps' ) && ~isempty( pbm.congrps ) )  % Check there are constraints!
         varargout{1} = evalHJv( 'Jv', pbm.congrps( varargin{4} ), varargin{2}, varargin{3}, [], pbm );
      else
         disp( ' ERROR: no constraint!' )
      end
   case 'cIJtxv'                                        % varargin = { pbm, x, v, I }
      if ( isfield( pbm, 'congrps' ) && ~isempty( pbm.congrps ) )  % Check there are constraints!
         varargout{1} = evalHJv( 'Jtv', pbm.congrps( varargin{4} ), varargin{2}, varargin{3}, [], pbm );
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The next function consider all problems in list_of_matlab_problems (whose files are in
%   the ./matlab_problems directory) and selects those whose SIF classification matches
%   that given by the input string varargin{1}. Matching is in the sense of regular expressions
%   (regexp).
%   If varargin{2} is empty (i.e. only varargin{1} is used as input argument), the function
%   prints the list of matching problems on standard output. Otherwise, the list is output in
%   the file whose name is a string passed as varargin{1} (Further components of varargin are
%   ignored).
%
%   If varargin{1} is 'help' or 'h', a message is printed on the standard output
%   describing the SIF classification scheme and an brief explanation of how to use the select
%   tool.
%
%   Thanks to Greta Malaspina (Firenze) for an inital implementation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'select'

   classif = varargin{1};
   if ( strcmp( classif, 'help') || strcmp( classif, 'h' ) )

      disp( '  ' )
      disp( ' === The classification scheme ===' )
      disp( '  ' )
      disp( ' A problem is classified by a string of the form' )
      disp( '    X-XXXXr-XX-n-m' )
      disp( ' The first character in the string identifies the problem collection' )
      disp( ' from which the problem is extracted. Possible values are' )
      disp( '    C the CUTEst collection;' )
      disp( '    S the SPARCO collection; and' )
      disp( '    N none of the above.' )
      disp( ' The character immediately following the first hyphen specifies the type' )
      disp( ' of variables occurring in the problem. Its possible values are' )
      disp( '    C the problem has continuous variables only;' )
      disp( '    I the problem has integer variables only;' )
      disp( '    B the problem has binary variables only; and' )
      disp( '    M the problem has variables of different types.' )
      disp( ' The second character after the first hyphen defines the type' )
      disp( ' of the problem''s objective function. Its possible values are' )
      disp( '    N no objective function is defined;' )
      disp( '    C the objective function is constant;' )
      disp( '    L the objective function is linear;' )
      disp( '    Q the objective function is quadratic;' )
      disp( '    S the objective function is a sum of squares; and' )
      disp( '    O the objective function is none of the above.' )
      disp( ' The third character after the first hyphen defines the type of constraints of the' )
      disp( ' problem. Its possible values are' )
      disp( '    U the problem is unconstrained;' )
      disp( '    X the problem’s only constraints are fixed variables;' )
      disp( '    B the problem’s only constraints are bounds on the variables;' )
      disp( '    N the problem’s constraints represent the adjacency matrix of a (linear)' )
      disp( '      network;' )
      disp( '    L the problem’s constraints are linear;' )
      disp( '    Q the problem’s constraints are quadratic; and' )
      disp( '    O the problem’s constraints are more general than any of the above alone.' )
      disp( ' The fourth character after the first hyphen indicates the smoothness of the problem.' )
      disp( ' There are two possible choices' )
      disp( '    R the problem is regular, that is, its first and second derivatives ' )
      disp( '      exist and are continuous everywhere; or' )
      disp( '    I the problem is irregular.' )
      disp( ' The integer (r) which corresponds to the fourth character of the string is' )
      disp( ' the degree of the highest derivatives provided analytically within the problem' )
      disp( ' description. It is restricted to being one of the single characters O, 1, or 2.' )
      disp( ' The character immediately following the second hyphen indicates the primary' )
      disp( ' origin and/or interest of the problem. Its possible values are' )
      disp( '    A the problem is academic, that is, has been constructed specifically by' )
      disp( '      researchers to test one or more algorithms;' )
      disp( '    M the problem is part of a modeling exercise where the actual value of the' )
      disp( '      solution is not used in a genuine practical application; and' )
      disp( '    R the problem’s solution is (or has been) actually used in a real')
      disp( '      application for purposes other than testing algorithms.' )
      disp( ' The next character in the string indicates whether or not the problem' )
      disp( ' description contains explicit internal variables. There are two possible' )
      disp( ' values, namely,' )
      disp( '    Y the problem description contains explicit internal variables; or' )
      disp( '    N the problem description does not contain any explicit internal variables.' )
      disp( ' The symbol(s) between the third and fourth hyphen indicate the number of' )
      disp( ' variables in the problem. Possible values are' )
      disp( '    V the number of variables in the problem can be chosen by the user; or' )
      disp( '    n a positive integer giving the actual (fixed) number of problem variables.' )
      disp( ' The symbol(s) after the fourth hyphen indicate the number of constraints' )
      disp( ' (other than fixed variables and bounds) in the problem. Note that fixed' )
      disp( ' variables are not considered as general constraints here. The two possible' )
      disp( ' values are' )
      disp( '    V the number of constraints in the problem can be chosen by the user; or' )
      disp( '    m a nonnegative integer giving the actual (fixed) number of constraints.' )
      disp( '  ' )
      disp( ' === Using the problem selection tool ===' )
      disp( '  ' )
      disp( ' To use the selection tool, you should first open Matlab in the parent of the' )
      disp( ' directory containing the Matlab problem files. The problem selection tool' )
      disp( ' may then by used in Matlab by calling s2mpjlib with action = ''select''. ' )
      disp( ' The second argument in such a call is a string specifying the class of problems' )
      disp( ' of interest.  This string is constructed by replacing by a dot each character' )
      disp( ' in the classification string for which all possible values are acceptable ' )
      disp( ' (the dot is a wildcard character).  For instance' )
      disp( '    s2mpjlib(''select'', ''C-CSU..-..-2-0'' ) ')
      disp( ' lists all CUTEst unconstrained sum-of-squares problems in two continuous' )
      disp( ' variables, while' )
      disp( '    s2mpjlib(''select'', ''C-C....-..-V-V'' ) ' )
      disp( ' lists all CUTEstproblems with variable number of continuous variables and' )
      disp( 'variable  number of constraints.' )
      disp( ' The classification strings ''unconstrained'', ''bound-constrained'', ' )
      disp( ' ''fixed-variables'', ''general-constraints'', ''variable-n'' and ' )
      disp( ' ''variable-m'' are also allowed.' )
      disp( ' NOTE: any regular expression may be used as the first argument of select ' )
      disp( '       to specify the problem class, so that, for instance, the previous ' )
      disp( '       selection can also be achieved by s2mpjlib(''select'', ''C-C.*V-V'' ) ')
      disp( ' Writing the list of selected problems to a file is obtained by specifying' )
      disp( ' the name of the file as a third argument of s2mplib, as in ')
      disp( '    s2mpjlib( ''select'', ''C-C....-..-V-V'', filename )' )
      
   else

      %  Modify the filter to cope with fixed numbers of variables/constraints with more
      %  than one digit.

      switch( classif )
      case 'unconstrained'
         classif = '.-..U.*';
      case 'bound-constrained'
         classif = '.-..B.*';
      case 'fixed-variables'
         classif = '.-..X.*';
      case 'general-constraints'
         classif = '.-..[LNQO].*';
      case 'variable-n'
         classif = '.-..B..-..-V-[V0-9]*';
      case 'variable-m'
         classif = '.-..B..-..-[V0-9]*-V';
      otherwise
         posh = strfind( classif, '-' );
         if ( length( posh ) >= 3 )
            oclassif = classif;
            if ( length( classif ) > posh(3) && classif( posh(3)+1 ) == '.' )
               classif = [ classif(1:posh(3)), '[V0-9]*' ];
               if ( length( oclassif ) > posh(3)+1 )
                  classif = [ classif, oclassif(posh(3)+2:end) ];
               end
            end
            if ( classif(end) == '.' )
               classif = [ classif(1:end-1), '[V0-9]*' ];
            end
         end
      end
      filter = [ 'classification = ''', classif ];

      %  Loop on the problems.

      list_of_problems = './list_of_matlab_problems';
      matlab_problems  = './matlab_problems/';

      if ( length( varargin ) > 1 )
         fid = fopen( varargin{2}, 'wt' );
      else
         fid = 1;
      end

      allprobs = readlines( list_of_problems );
      for i=1:length( allprobs )
         theprob = allprobs(i);
         problem = matlab_problems + theprob;
         if ( isfile( problem ) )
            if ( ~isempty( regexp( fileread( problem ), filter, 'match' ) ) )
               fprintf( fid, '%s\n', theprob );
            end
         end
      end
      if ( fid ~= 1 )
         fclose( fid );
      end
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
   m    = length( glist );
   cx   = sparse( m, 1 );
   ic   = 0;
   if ( isfield( pbm, 'conderlvl' ) )
      lder = length( pbm.conderlvl );
   end
end
if ( nargout > 1 )
   if ( isobj )
      if ( ~isfield( pbm, 'objderlvl' ) || pbm.objderlvl >= 1 )
         if ( redprec )
            gx = 0*x;
         else
            gx = sparse( n,1 );
         end
      else
         gx = sparse( n,1 );
         gx( 1 ) = NaN;
      end
   else
       if ( ~isfield( pbm, 'conderlvl' ) || any ( pbm.conderlvl(1) >= 1 ) )
         if ( redprec )
            Jx = sparse( m, n )*(0*x(1));
         else
            Jx = sparse( m, n );
         end
      else
         Jx = sparse( m, n );
         Jx( 1, 1 ) =  NaN;
      end
   end
   if ( nargout > 2 )
      if ( isobj )
         if ( ~isfield( pbm, 'objderlvl' ) || pbm.objderlvl >= 2 )
            if ( redprec )
               Hx = sparse( n, n )*(0*x(1));
            else
               Hx = sparse( n, n );
            end
         else
            Hx = sparse( n, n );
            Hx( 1, 1 ) = NaN;
         end
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

   %  Find the level of available derivatives for the group.

   if (isobj )
      if ( ~isfield( pbm, 'objderlvl' ) )
         derlvl = 2;
      else
         derlvl = pbm.objderlvl;
      end
   else
      if ( isfield( pbm, 'conderlvl' ) )
         if ( lder == 1 )
            derlvl = pbm.conderlvl( 1 );
         else
            [ ~, posg ] = ismember( ig, pbm.congrps );
            derlvl      = pbm.conderlvl( posg );
         end
      else
         derlvl = 2;
      end
   end
   nout = min( nargout, derlvl + 1 );

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
         gin = sparse( n, 1 ) *(0*x(1));
      else
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
         switch ( nout )

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
            if ( derlvl >= 1 )
               gx = gx + grada * gin /gsc ;
            else
               gx = sparse( n, 1 );
               gx(1) = NaN;
            end
         case 3
           [ fa, grada, Hessa ] = feval( pbm.name, egname, fin, ig );
            fx   = fx + fa / gsc;
            if ( derlvl >= 1 )
               gx   = gx + grada * gin / gsc;
            else
               gx = Nan*ones( n, 1 );
            end
            if ( derlvl >= 2 )            
               Hx = Hx + ( ( Hessa * gin  ) * gin.' + grada * Hin ) / gsc;
            else
               Hx = pqarse( n, n );
               Hx( 1, 1 ) = NaN;
            end
         end
      else
         ic = ic + 1;
         switch ( nargout )
         case 1
            cx( ic ) = feval( pbm.name, egname, fin, ig ) / gsc;
         case 2
            [ fa, grada ] = feval( pbm.name, egname, fin, ig );
            cx( ic )      = fa / gsc;
            if ( derlvl >= 1 )
               Jx(ic,1:n) = grada * gin.' / gsc;
            else
               Jx(ic,1:n) = deal( NaN );
            end            
         case 3
            [ fa, grada, Hessa ] = feval( pbm.name, egname, fin, ig );
            cx( ic )   = fa / gsc;
            if ( derlvl >= 1 )
               Jx(ic,1:n) = grada * gin.' / gsc;
            else
               Jx(ic,1:n) = deal( NaN );
            end            
            if ( derlvl >= 2 )
               Hx{end+1} = ( ( Hessa * gin ) * gin.' + grada * Hin ) / gsc;
            else
               Hx{end+1} = sparse( n, n );
               Hx{end}( 1, 1 ) = NaN;
            end
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
            if ( derlvl >= 1 )
               gx = gx + gin / gsc;
            else
               gx = sparse( n, 1 );
               gx( 1 ) = NaN;
            end
         case 3
            fx = fx + fin / gsc;
            if ( derlvl >= 1 )
               gx = gx + gin / gsc;
            else
               gx = NaN * sparse( n, 1 );
            end
            if ( derlvl >= 2 )
                Hx = Hx + Hin / gsc;
            else
                Hx = sparse( n, n );
                Hx( 1, 1 ) = NaN;
            end
         end
      else
         ic = ic + 1;
         switch ( nargout )
         case 1
            cx( ic )   = fin  / gsc;
         case 2
            cx( ic )   = fin  / gsc;
            if ( derlvl >= 1 )
               Jx(ic,1:n) = gin.'/ gsc;
            else
               Jx(ic,1:n) = deal( NaN);
            end
         case 3
            cx( ic )   = fin  / gsc;
            if ( derlvl >= 1 )
               Jx(ic,1:n) = gin.'/ gsc;
            else
               Jx(ic,1:n) = deal( NaN);
            end
            if ( derlvl >= 2 )
               Hx{end+1}  = Hin  / gsc;
            else
               Hx{end+1}  = sparse( n, n );
               Hx{end}( 1, 1 ) = NaN;
            end
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
%   mode = 'Hv'  : evaluate the product of the objective's Hessian times v (glist = obj groups)
%   mode = 'HIv' : evaluate the product of the constraints' Hessian  times v time the multiplier y
%                 (glist = cons groups)
%   mode = 'Jv'  : evaluate the product of the constraints' Jacobian times v (glist = cons groups)
%   The vector y is unused (and unreferenced) for modes 'Hv' and 'Jv'.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

debug = 0;%D

%  Avoid computing anything when the relevant derivatives are missing.

n = length( x );
switch ( mode )
case 'Hv'
   if ( isfield( pbm, 'objcondlvl' ) && pbm.objderlvl < 2 )
      HJv = sparse( n, 1 );
      HJv( 1, 1 ) = NaN;
      return
   else
      HJv    = zeros( n, 1 );
      derlvl = 2;
   end
case 'HIv'
   m      = length( glist );
   lder   = length( pbm.conderlvl );
   if ( isfield( pbm, 'conderlvl' ) && any( pbm.conderlvl < 2 ) )
      HJv = sparse( n, 1 );
      HJv( 1, 1 ) = NaN;
      return
   else
      HJv    = zeros( n, 1 );
      ic     = 0;
      derlvl = 2;
   end
case { 'Jv', 'Jtv' }
   m      = length( glist );
   lder   = length( pbm.conderlvl );
   if ( isfield( pbm, 'conderlvl' ) && any( pbm.conderlvl < 1 ) )
      switch ( mode )
      case 'Jv'
         HJv = sparse( m, 1 );
      case 'Jtv'
         HJv = sparse( n, 1 );
      end
      HJv( 1, 1 ) = NaN;
      return
   else
      switch ( mode )
      case 'Jv'
         HJv = zeros( m, 1 );
      case 'Jtv'
         HJv = zeros( n, 1 );
      end
      ic  = 0;
   end
otherwise
   disp( [' ERROR: unknown mode', mode, ' at entry of evalHJv.' ] )
end

%  See if reduced precision is requested.

if ( isfield( pbm, 'ndigs' ) )
   ndigs   = pbm.ndigs;
   redprec = 1;
else
   redprec = 0;
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

   %  Find the level of available derivatives for the group.

   if ( ismember( mode, { 'Jv', 'Jtv' } ) )
      if ( isfield( pbm, 'conderlvl' ) )
         if ( lder == 1 )
            derlvl = pbm.conderlvl( 1 );
         else
            [ ~, posg ] = ismember( ig, pbm.congrps );
            derlvl      = pbm.conderlvl( posg );
         end
      else
         derlvl = 2;
      end
      
      %  Avoid computation for group ig if its first derivative is missing.
      
      if ( derlvl < 1 )
         ic      = ic + 1;
         HJv(ic) = NaN;
         continue
      end
   end

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
   
   if( debug )
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
         
         elseif ( derlvl >= 1 )
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
      if ( derlvl >= 1 )
         switch ( mode )
         case 'Jv'
            if ( strcmp( egname, 'TRIVIAL' ) )
               HJv(ic) = gin.' * v / gsc;
            else
               [ ~, grada ] = feval( pbm.name, egname, fin, ig );
               HJv(ic)      = grada * gin.' * v / gsc;
            end
         case 'Jtv'
            if ( strcmp( egname, 'TRIVIAL' ) )
               HJv = HJv + gin * v( ic ) / gsc;
            else
               [ ~, grada ] = feval( pbm.name, egname, fin, ig );
               HJv = HJv + grada * gin * v( ic ) / gsc;
            end
         end
      else
         HJv(ic) = NaN;
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
%   Hessian times a vector v, depending on'action'.
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

if( ~isempty( gconlist ) )
   switch( nargout )
   case 1
      switch( action )
      case { 'Lxy', 'Lgxy', 'LgHxy', 'LIxy', 'LIgxy', 'LIgHxy' }
         varargout{1}    = varargout{1} + y' * evalgrsum( 0, gconlist, x, pbm );
      case { 'LHxyv', 'LIHxyv' }
         varargout{1} = varargout{1} + evalHJv( 'HIv', gconlist, x, v, y, pbm );
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
end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
