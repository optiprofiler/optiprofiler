function varargout = SIMPLEU( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : SIMPLEU
%    --------
%
%    A S2X/Matlab template allowing users to input their (unstructured)
%    unconstrained problem easily.
%
%    It is assumed thet the user has available a Matlab function called
%    your_obj.function.m defining an objective function to minimize, whose
%    interface is
%
%          [ varargout{1:nargout} ] = your_obj_function( x, nargout )
%
%    where
%           x        is a REAL vector of variable's values
%           nargout  is an integer such that
%                    if nargout = 1, then varargout{1} is the value of
%                                    the function evaluated at x
%                    if nargout = 2; then, in addition, varargout{2} is
%                                    the value of the function's gradient
%                                    evaluated at x
%                    if nargout = 3; then, in addition, varargout{3} is
%                                    the value of the function's Hessian
%                                    evaluated at x (as a sparse matrix).
%
%    The template then provides the interface between the user's objective
%    function and the S2X environment, allowing the use of the s2xlib.m
%    library for ecvaluations.
% 
%    classification = 'OUR2-AN-V-0'
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'SIMPLEU';

switch(action)

    case 'setup'
        %%%%%%%%%%%%%%%%% NUMBER OF VARIABLES %%%%%%%%%%%%%%
        pb.n = varargin{1};
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% BOUNDS %%%%%%%%%%%%%%%%%%%%%%%
        pb.xlower      = -Inf*ones(pb.n,1);
        pb.xupper      = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%% DON'T MODIFY  %%%%%%%%%%%%%%%%%%%
        pb.name        = name;                             %
        pb.sifpbname   = name;                             %
        pbm.name       = name;                             %
        pbm.objgrps    = [1];                              %
        pb.m           = 0;                                %
        pbm.elftype{1} = 'ETYPE';                          %
        pbm.elvar{1}   = [1:pb.n];                         %
        pbm.grelt{1}   = [1];                              %
        pbm.gconst     = 0.0;                              %
        pb.pbclass     = 'OUR2-AN-V-0';                    %
        pb.xnames = {};                                    %
        for i = 1:pb.n                                     %
            pb.xnames{i} = ['X',int2str(i)];               %
        end                                                %
        varargout{1}   = pb;                               %
        varargout{2}   = pbm;                              %

    %%%%%%%%%%%%%%%% THE NONLINEAR OBJECTIVE %%%%%%%%%%%%%%%

    case 'ETYPE'                                           %

        x  = varargin{1};                                  %
        [varargout{1:nargout}] = your_obj_function( x, nargout );   

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy','LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [0,0];
            [varargout{1:max(1,nargout)}] = s2xlib(action,pbm,varargin{:});
        else
            disp(['ERROR: please run ',name,' with action = setup'])
            [varargout{1:nargout}] = deal(repmat(NaN,1:nargout));
            end

    otherwise
        disp([' ERROR: unknown action ',action,' requested from ',name,'.m'])
    end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

