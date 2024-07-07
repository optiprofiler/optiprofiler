function varargout = LUKSAN17(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : LUKSAN17
%    *********
% 
%    Problem 17 (sparse trigonometric) in the paper
% 
%      L. Luksan
%      Hybrid methods in large sparse nonlinear least squares
%      J. Optimization Theory & Applications 89(3) 575-595 (1996)
% 
%    SIF input: Nick Gould, June 2017.
% 
%    classification = 'NOR2-AN-V-V'
% 
%   seed for dimensions
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'LUKSAN17';

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
        v_('S') = 49;
        v_('N') = 2*v_('S');
        v_('N') = 2+v_('N');
        v_('M') = 4*v_('S');
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('4') = 4;
        v_('ONE') = 1.0;
        v_('Y1') = 30.6;
        v_('Y2') = 72.2;
        v_('Y3') = 124.4;
        v_('Y4') = 187.4;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('M')
            [ig,ig_] = s2mpjlib('ii',['E',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['E',int2str(I)];
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
        legrps = find(strcmp(gtype,'<='));
        eqgrps = find(strcmp(gtype,'=='));
        gegrps = find(strcmp(gtype,'>='));
        pb.nle = length(legrps);
        pb.neq = length(eqgrps);
        pb.nge = length(gegrps);
        pb.m   = pb.nle+pb.neq+pb.nge;
        pbm.congrps = [ legrps, eqgrps, gegrps ];
        [pb.cnames{1:pb.m}] = deal(cnames{pbm.congrps});
        pb.nob = ngrp-pb.m;
        pbm.objgrps = find(strcmp(gtype,'<>'));
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        v_('K') = 1;
        for J=v_('1'):v_('S')
            pbm.gconst(ig_(['E',int2str(round(v_('K')))])) = v_('Y1');
            v_('K') = 1+v_('K');
            pbm.gconst(ig_(['E',int2str(round(v_('K')))])) = v_('Y2');
            v_('K') = 1+v_('K');
            pbm.gconst(ig_(['E',int2str(round(v_('K')))])) = v_('Y3');
            v_('K') = 1+v_('K');
            pbm.gconst(ig_(['E',int2str(round(v_('K')))])) = v_('Y4');
            v_('K') = 1+v_('K');
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        for I=v_('1'):v_('4'):v_('N')
            pb.x0(ix_(['X',int2str(I)]),1) = -0.8;
        end
        for I=v_('2'):v_('4'):v_('N')
            pb.x0(ix_(['X',int2str(I)]),1) = 1.2;
        end
        for I=v_('3'):v_('4'):v_('N')
            pb.x0(ix_(['X',int2str(I)]),1) = -1.2;
        end
        for I=v_('4'):v_('4'):v_('N')
            pb.x0(ix_(['X',int2str(I)]),1) = 0.8;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eACOSX',iet_);
        elftv{it}{1} = 'X';
        elftp{it}{1} = 'A';
        [it,iet_] = s2mpjlib( 'ii', 'eASINX',iet_);
        elftv{it}{1} = 'X';
        elftp{it}{1} = 'A';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for Q=v_('1'):v_('4')
            v_('RQ') = Q;
            v_('RQ2') = v_('RQ')*v_('RQ');
            v_('K') = 1;
            v_('I') = 0;
            for J=v_('1'):v_('S')
                v_('I+Q') = v_('I')+Q;
                for L=v_('1'):v_('4')
                    v_('RL') = L;
                    v_('RL2') = v_('RL')*v_('RL');
                    v_('A') = v_('RL')*v_('RQ2');
                    v_('A') = -1.0*v_('A');
                    ename = ['S',int2str(round(v_('K'))),',',int2str(Q)];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'eASINX';
                    ielftype(ie) = iet_('eASINX');
                    ename = ['S',int2str(round(v_('K'))),',',int2str(Q)];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    vname = ['X',int2str(round(v_('I+Q')))];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('X',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    ename = ['S',int2str(round(v_('K'))),',',int2str(Q)];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    [~,posep] = ismember('A',elftp{ielftype(ie)});
                    pbm.elpar{ie}(posep) = v_('A');
                    v_('A') = v_('RL2')*v_('RQ');
                    ename = ['C',int2str(round(v_('K'))),',',int2str(Q)];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    pbm.elftype{ie} = 'eACOSX';
                    ielftype(ie) = iet_('eACOSX');
                    ename = ['C',int2str(round(v_('K'))),',',int2str(Q)];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    vname = ['X',int2str(round(v_('I+Q')))];
                    [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
                    posev = find(strcmp('X',elftv{ielftype(ie)}));
                    pbm.elvar{ie}(posev) = iv;
                    ename = ['C',int2str(round(v_('K'))),',',int2str(Q)];
                    [ie,ie_] = s2mpjlib('ii',ename,ie_);
                    [~,posep] = ismember('A',elftp{ielftype(ie)});
                    pbm.elpar{ie}(posep) = v_('A');
                    v_('K') = 1+v_('K');
                end
                v_('I') = 2+v_('I');
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for K=v_('1'):v_('M')
            for Q=v_('1'):v_('4')
                ig = ig_(['E',int2str(K)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['S',int2str(K),',',int2str(Q)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = 1.;
                posel = posel+1;
                pbm.grelt{ig}(posel) = ie_(['C',int2str(K),',',int2str(Q)]);
                pbm.grelw{ig}(posel) = 1.;
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN                0.0
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'NOR2-AN-V-V';
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end
% **********************
%  SET UP THE FUNCTION *
%  AND RANGE ROUTINES  *
% **********************

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eASINX'

        EV_  = varargin{1};
        iel_ = varargin{2};
        ASINX = pbm.elpar{iel_}(1)*sin(EV_(1));
        varargout{1} = ASINX;
        if(nargout>1)
            g_(1,1) = pbm.elpar{iel_}(1)*cos(EV_(1));
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = -ASINX;
                varargout{3} = H_;
            end
        end

    case 'eACOSX'

        EV_  = varargin{1};
        iel_ = varargin{2};
        ACOSX = pbm.elpar{iel_}(1)*cos(EV_(1));
        varargout{1} = ACOSX;
        if(nargout>1)
            g_(1,1) = -pbm.elpar{iel_}(1)*sin(EV_(1));
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = -ACOSX;
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

