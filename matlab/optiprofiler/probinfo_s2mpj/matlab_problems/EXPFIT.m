function varargout = EXPFIT(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : EXPFIT
%    *********
%    A simple exponential fit in 2 variables
% 
%    Source:
%    A.R. Conn, N. Gould and Ph.L. Toint,
%    "LANCELOT, a Fortran package for large-scale nonlinear optimization",
%    Springer Verlag, FUNDP, 1992.
% 
%    SIF input: Ph. Toint, Jan 1991.
% 
%    classification = 'SUR2-AN-2-0'
% 
%    Number of points
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'EXPFIT';

switch(action)

    case {'setup','setup_redprec'}

        if(isfield(pbm,'ndigs'))
            rmfield(pbm,'ndigs');
        end
        if(strcmp(action,'setup_redprec'))
            pbm.ndigs = max(1,min(15,varargin{end}));
            nargs     = nargin-2;
        else
            nargs = nargin-1;
        end
        pb.name   = name;
        pbm.name  = name;
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        v_('P') = 10;
        v_('H') = 0.25;
        v_('1') = 1;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','ALPHA',ix_);
        pb.xnames{iv} = 'ALPHA';
        [iv,ix_] = s2mpjlib('ii','BETA',ix_);
        pb.xnames{iv} = 'BETA';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for i=v_('1'):v_('P')
            [ig,ig_] = s2mpjlib('ii',['R',int2str(i)],ig_);
            gtype{ig} = '<>';
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for i=v_('1'):v_('P')
            v_('Reali') = i;
            v_('Reali*H') = v_('Reali')*v_('H');
            pbm.gconst(ig_(['R',int2str(i)])) = v_('Reali*H');
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eEXPIH',iet_);
        elftv{it}{1} = 'V';
        elftv{it}{2} = 'W';
        elftp{it}{1} = 'RI';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for i=v_('1'):v_('P')
            v_('Reali') = i;
            ename = ['E',int2str(i)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eEXPIH';
            ielftype(ie) = iet_('eEXPIH');
            vname = 'ALPHA';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'BETA';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('RI',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('Reali');
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for ig = 1:ngrp
            pbm.grftype{ig} = 'gL2';
        end
        for i=v_('1'):v_('P')
            ig = ig_(['R',int2str(i)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(i)]);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'SUR2-AN-2-0';
        pb.x0          = zeros(pb.n,1);
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eEXPIH'

        EV_  = varargin{1};
        iel_ = varargin{2};
        IH = 0.25*pbm.elpar{iel_}(1);
        EXPWIH = exp(EV_(2)*IH);
        varargout{1} = EV_(1)*EXPWIH;
        if(nargout>1)
            g_(1,1) = EXPWIH;
            g_(2,1) = EV_(1)*IH*EXPWIH;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = IH*EXPWIH;
                H_(2,1) = H_(1,2);
                H_(2,2) = EV_(1)*IH*IH*EXPWIH;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gL2'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_^2;
        if(nargout>1)
            g_ = 2.0*GVAR_;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 2.0;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

