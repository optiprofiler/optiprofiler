function varargout = ROSSIMP1NE(action,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Problem : ROSSIMP1NE
%    The ever famous 2 variables Rosenbrock "banana valley" problem
%    This version considers the problem as a system of nonlinear equations.
%    classification = 'NOR2-AN-2-0'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
persistent pbm;
name = 'ROSSIMP1NE';
switch(action)
    case 'setup'
        pb.name      = name;
        pbm.name     = name;
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(2,2);
        pbm.A(1,2) = 1.0;
        pbm.gscale(1) = 0.1;
        pbm.A(2,1) = 1.0;
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = 2;
        pbm.congrps = [ 1, 2 ];
        pb.neq      = 2;
        pb.m        = 2;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst =[ 0.0; 1 ];
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = [ -Inf; -Inf ];
        pb.xupper = [ +Inf; +Inf ];
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = [ -1.2; 1.0 ];
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        pbm.elftype{1} = 'eSQ';
        pbm.elvar{1}   = [ 1 ];
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        pbm.grelt{1} = [ 1 ];
        pbm.grelw{1} = [ -1.0 ];
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'NOR2-AN-2-0';
        varargout{1} = pb;
        varargout{2} = pbm;
    case 'eSQ'
        EV_  = varargin{1};
        varargout{1} = EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = EV_(1)+EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0;
                varargout{3} = H_;
            end
        end
    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%
    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy','LHxyv','LIHxyv'}
        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [0,0];
            [varargout{1:max(1,nargout)}] = s2mpjlib(action,pbm,varargin{:});
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

