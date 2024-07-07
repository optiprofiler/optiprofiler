function varargout = LUKSAN14(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : LUKSAN14
%    *********
% 
%    Problem 14 (chained and modified HS53) in the paper
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

name = 'LUKSAN14';

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
        v_('S') = 32;
        v_('N') = 3*v_('S');
        v_('N') = 2+v_('N');
        v_('M') = 7*v_('S');
        v_('1') = 1;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        v_('I') = 1;
        v_('K') = 1;
        for J=v_('1'):v_('S')
            v_('K+1') = 1+v_('K');
            v_('K+2') = 2+v_('K');
            v_('K+3') = 3+v_('K');
            v_('K+4') = 4+v_('K');
            v_('K+5') = 5+v_('K');
            v_('K+6') = 6+v_('K');
            v_('I+1') = 1+v_('I');
            v_('I+2') = 2+v_('I');
            v_('I+3') = 3+v_('I');
            v_('I+4') = 4+v_('I');
            [ig,ig_] = s2mpjlib('ii',['E',int2str(round(v_('K')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['E',int2str(round(v_('K')))];
            iv = ix_(['X',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -10.0e0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -10.0e0;
            end
            [ig,ig_] = s2mpjlib('ii',['E',int2str(round(v_('K+1')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['E',int2str(round(v_('K+1')))];
            iv = ix_(['X',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0e0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0e0;
            end
            iv = ix_(['X',int2str(round(v_('I+2')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0e0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0e0;
            end
            [ig,ig_] = s2mpjlib('ii',['E',int2str(round(v_('K+2')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['E',int2str(round(v_('K+2')))];
            iv = ix_(['X',int2str(round(v_('I+3')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0e0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0e0;
            end
            [ig,ig_] = s2mpjlib('ii',['E',int2str(round(v_('K+3')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['E',int2str(round(v_('K+3')))];
            iv = ix_(['X',int2str(round(v_('I+4')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0e0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0e0;
            end
            [ig,ig_] = s2mpjlib('ii',['E',int2str(round(v_('K+4')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['E',int2str(round(v_('K+4')))];
            iv = ix_(['X',int2str(round(v_('I')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0e0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0e0;
            end
            iv = ix_(['X',int2str(round(v_('I+1')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 3.0e0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 3.0e0;
            end
            [ig,ig_] = s2mpjlib('ii',['E',int2str(round(v_('K+5')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['E',int2str(round(v_('K+5')))];
            iv = ix_(['X',int2str(round(v_('I+2')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0e0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0e0;
            end
            iv = ix_(['X',int2str(round(v_('I+3')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0e0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0e0;
            end
            [ig,ig_] = s2mpjlib('ii',['E',int2str(round(v_('K+5')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['E',int2str(round(v_('K+5')))];
            iv = ix_(['X',int2str(round(v_('I+4')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -2.0e0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -2.0e0;
            end
            [ig,ig_] = s2mpjlib('ii',['E',int2str(round(v_('K+6')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['E',int2str(round(v_('K+6')))];
            iv = ix_(['X',int2str(round(v_('I+4')))]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -10.0e0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -10.0e0;
            end
            v_('I') = 3+v_('I');
            v_('K') = 7+v_('K');
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
            v_('K+1') = 1+v_('K');
            v_('K+2') = 2+v_('K');
            v_('K+3') = 3+v_('K');
            pbm.gconst(ig_(['E',int2str(round(v_('K+1')))])) = 2.0;
            pbm.gconst(ig_(['E',int2str(round(v_('K+2')))])) = 1.0;
            pbm.gconst(ig_(['E',int2str(round(v_('K+3')))])) = 1.0;
            v_('K') = 7+v_('K');
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        for I=v_('1'):v_('N')
            pb.x0(ix_(['X',int2str(I)]),1) = -1.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQR',iet_);
        elftv{it}{1} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        v_('I') = 1;
        v_('K') = 1;
        for J=v_('1'):v_('S')
            v_('K+6') = 6+v_('K');
            v_('I+1') = 1+v_('I');
            ename = ['E',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
            ename = ['E',int2str(round(v_('K')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('I')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['E',int2str(round(v_('K+6')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
            ename = ['E',int2str(round(v_('K+6')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('I+1')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            v_('I') = 3+v_('I');
            v_('K') = 7+v_('K');
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        v_('K') = 1;
        for J=v_('1'):v_('S')
            v_('K+6') = 6+v_('K');
            ig = ig_(['E',int2str(round(v_('K')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('K')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 10.0;
            ig = ig_(['E',int2str(round(v_('K+6')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('K+6')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 10.0;
            v_('K') = 7+v_('K');
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

    case 'eSQR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(1);
        if(nargout>1)
            g_(1,1) = EV_(1)+EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = 2.0e0;
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

