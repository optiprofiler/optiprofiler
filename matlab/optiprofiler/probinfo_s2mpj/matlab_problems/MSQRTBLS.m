function varargout = MSQRTBLS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : MSQRTBLS
%    *********
% 
%    The dense matrix square root problem by Nocedal and Liu (Case 1)
% 
%    This is a least-squares variant of problem MSQRTB.
% 
%    Source:  problem 204 (p. 93) in
%    A.R. Buckley,
%    "Test functions for unconstrained minimization",
%    TR 1989CS-3, Mathematics, statistics and computing centre,
%    Dalhousie University, Halifax (CDN), 1989.
% 
%    SIF input: Ph. Toint, Dec 1989.
% 
%    classification = 'SUR2-AN-V-V'
% 
%    Dimension of the matrix ( at least 3)
% 
%       Alternative values for the SIF file parameters:
% IE P                   3              $-PARAMETER n = 9     original value
% IE P                   7              $-PARAMETER n = 49
% IE P                   10             $-PARAMETER n = 100
% IE P                   23             $-PARAMETER n = 529
% IE P                   32             $-PARAMETER n = 1024
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'MSQRTBLS';

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
        if(nargs<1)
            v_('P') = 5;  %  SIF file default value
        else
            v_('P') = varargin{1};
        end
% IE P                   70             $-PARAMETER n = 4900
        v_('N') = v_('P')*v_('P');
        v_('1') = 1;
        v_('K') = 0.0;
        for I=v_('1'):v_('P')
            for J=v_('1'):v_('P')
                v_('K') = 1.0+v_('K');
                v_('K2') = v_('K')*v_('K');
                v_(['B',int2str(I),',',int2str(J)]) = sin(v_('K2'));
            end
        end
        v_('B3,1') = 0.0;
        for I=v_('1'):v_('P')
            for J=v_('1'):v_('P')
                v_(['A',int2str(I),',',int2str(J)]) = 0.0;
                for T=v_('1'):v_('P')
                    v_('PROD') =...
                          v_(['B',int2str(I),',',int2str(T)])*v_(['B',int2str(T),',',int2str(J)]);
                    v_(['A',int2str(I),',',int2str(J)]) = v_(['A',int2str(I),',',int2str(J)])+...
                         v_('PROD');
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('P')
            for J=v_('1'):v_('P')
                [iv,ix_] = s2mpjlib('ii',['X',int2str(I),',',int2str(J)],ix_);
                pb.xnames{iv} = ['X',int2str(I),',',int2str(J)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('P')
            for J=v_('1'):v_('P')
                [ig,ig_] = s2mpjlib('ii',['G',int2str(I),',',int2str(J)],ig_);
                gtype{ig} = '<>';
            end
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('1'):v_('P')
            for J=v_('1'):v_('P')
                pbm.gconst(ig_(['G',int2str(I),',',int2str(J)])) =...
                      v_(['A',int2str(I),',',int2str(J)]);
            end
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        v_('K') = 0.0;
        for I=v_('1'):v_('P')
            for J=v_('1'):v_('P')
                v_('K') = 1.0+v_('K');
                v_('K2') = v_('K')*v_('K');
                v_('SK2') = sin(v_('K2'));
                v_('-4SK2/5') = -0.8*v_('SK2');
                v_('XIJ') = v_(['B',int2str(I),',',int2str(J)])+v_('-4SK2/5');
                pb.x0(ix_(['X',int2str(I),',',int2str(J)]),1) = v_('XIJ');
            end
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'en2PR',iet_);
        elftv{it}{1} = 'XIT';
        elftv{it}{2} = 'XTJ';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('P')
            for J=v_('1'):v_('P')
                for T=v_('1'):v_('P')
                    ename = ['E',int2str(I),',',int2str(J),',',int2str(T)];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'en2PR';
                    ielftype(ie) = iet_('en2PR');
                    vname = ['X',int2str(I),',',int2str(T)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('XIT',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    vname = ['X',int2str(T),',',int2str(J)];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('XTJ',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                end
            end
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
        for I=v_('1'):v_('P')
            for J=v_('1'):v_('P')
                for T=v_('1'):v_('P')
                    ig = ig_(['G',int2str(I),',',int2str(J)]);
                    posel = length(pbm.grelt{ig})+1;
                    pbm.grelt{ig}(posel) = ie_(['E',int2str(I),',',int2str(J),',',int2str(T)]);
                    pbm.grelw{ig}(posel) = 1.;
                end
            end
        end
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'SUR2-AN-V-V';
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'en2PR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(2);
        if(nargout>1)
            g_(1,1) = EV_(2);
            g_(2,1) = EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = 1.0;
                H_(2,1) = H_(1,2);
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%%%% NONLINEAR GROUPS  %%%%%%%%%%%%%%%

    case 'gL2'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_*GVAR_;
        if(nargout>1)
            g_ = GVAR_+GVAR_;
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

