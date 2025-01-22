function varargout = LUKVLI11(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : LUKVLI11
%    *********
% 
%    Source: Problem 5.11, the chained HS46 problem, 
%    due to L. Luksan and J. Vlcek,
%    "Sparse and partially separable test problems for 
%    unconstrained and equality constrained optimization",
%    Technical Report 767, Inst. Computer Science, Academy of Sciences
%    of the Czech Republic, 182 07 Prague, Czech Republic, 1999
% 
%    Equality constraints changed to inequalities
% 
%    SIF input: Nick Gould, April 2001
% 
%    classification = 'C-COOR2-AY-V-V'
% 
%    some useful parameters, including N, the number of variables.
% 
%       Alternative values for the SIF file parameters:
% IE N                   98             $-PARAMETER
% IE N                   998            $-PARAMETER
% IE N                   9998           $-PARAMETER
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'LUKVLI11';

switch(action)

    case {'setup','setup_redprec'}

        if(strcmp(action,'setup_redprec'))
            if(isfield(pbm,'ndigs'))
                rmfield(pbm,'ndigs');
            end
            pbm.ndigs = max(1,min(15,varargin{end}));
            nargs     = nargin-2;
        else
            nargs = nargin-1;
        end
        pb.name   = name;
        pbm.name  = name;
        %%%%%%%%%%%%%%%%%%%%  PREAMBLE %%%%%%%%%%%%%%%%%%%%
        v_  = containers.Map('KeyType','char', 'ValueType', 'double');
        ix_ = containers.Map('KeyType','char', 'ValueType', 'double');
        ig_ = containers.Map('KeyType','char', 'ValueType', 'double');
        if(nargs<1)
            v_('N') = 8;  %  SIF file default value
        else
            v_('N') = varargin{1};
        end
% IE N                   99998          $-PARAMETER
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('4') = 4;
        v_('5') = 5;
        v_('N-2') = -2+v_('N');
        v_('(N-2)/3') = fix(v_('N-2')/v_('3'));
        v_('NC') = v_('2')*v_('(N-2)/3');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for I=v_('1'):v_('(N-2)/3')
            v_('I-1') = -1+I;
            v_('J') = v_('3')*v_('I-1');
            v_('J+1') = 1+v_('J');
            v_('J+2') = 2+v_('J');
            v_('J+3') = 3+v_('J');
            v_('J+4') = 4+v_('J');
            v_('J+5') = 5+v_('J');
            [ig,ig_] = s2mpjlib('ii',['OBJ1',int2str(I)],ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('J+1')))]);
            valA(end+1) = 1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('J+2')))]);
            valA(end+1) = -1.0;
            [ig,ig_] = s2mpjlib('ii',['OBJ2',int2str(I)],ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('J+3')))]);
            valA(end+1) = 1.0;
            [ig,ig_] = s2mpjlib('ii',['OBJ3',int2str(I)],ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('J+4')))]);
            valA(end+1) = 1.0;
            [ig,ig_] = s2mpjlib('ii',['OBJ4',int2str(I)],ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('J+5')))]);
            valA(end+1) = 1.0;
        end
        for K=v_('1'):v_('2'):v_('NC')
            v_('K+1') = 1+K;
            [ig,ig_] = s2mpjlib('ii',['C',int2str(K)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['C',int2str(K)];
            [ig,ig_] = s2mpjlib('ii',['C',int2str(round(v_('K+1')))],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['C',int2str(round(v_('K+1')))];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('K+1')))]);
            valA(end+1) = 1.0;
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
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
        for I=v_('1'):v_('(N-2)/3')
            pbm.gconst(ig_(['OBJ2',int2str(I)])) = 1.0;
            pbm.gconst(ig_(['OBJ3',int2str(I)])) = 1.0;
            pbm.gconst(ig_(['OBJ4',int2str(I)])) = 1.0;
        end
        for K=v_('1'):v_('2'):v_('NC')
            v_('K+1') = 1+K;
            pbm.gconst(ig_(['C',int2str(K)])) = 1.0;
            pbm.gconst(ig_(['C',int2str(round(v_('K+1')))])) = 2.0;
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        for I=v_('1'):v_('3'):v_('N')
            pb.x0(ix_(['X',int2str(I)]),1) = 2.0;
        end
        for I=v_('2'):v_('3'):v_('N')
            pb.x0(ix_(['X',int2str(I)]),1) = 1.5;
        end
        for I=v_('3'):v_('3'):v_('N')
            pb.x0(ix_(['X',int2str(I)]),1) = 0.5;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eC21',iet_);
        elftv{it}{1} = 'V';
        elftv{it}{2} = 'W';
        [it,iet_] = s2mpjlib( 'ii', 'eS',iet_);
        elftv{it}{1} = 'V';
        elftv{it}{2} = 'W';
        [it,iet_] = s2mpjlib( 'ii', 'eC42',iet_);
        elftv{it}{1} = 'V';
        elftv{it}{2} = 'W';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for K=v_('1'):v_('2'):v_('NC')
            v_('K+1') = 1+K;
            v_('K+2') = 2+K;
            v_('K+3') = 3+K;
            v_('K+4') = 4+K;
            ename = ['EA',int2str(K)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eC21';
            ielftype(ie) = iet_('eC21');
            vname = ['X',int2str(K)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('K+3')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['EB',int2str(K)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eS';
            ielftype(ie) = iet_('eS');
            vname = ['X',int2str(round(v_('K+3')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['X',int2str(round(v_('K+4')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['E',int2str(round(v_('K+1')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eC21';
            ielftype(ie) = iet_('eC21');
            ename = ['E',int2str(round(v_('K+1')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('K+2')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            ename = ['E',int2str(round(v_('K+1')))];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            vname = ['X',int2str(round(v_('K+3')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('W',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        [it,igt_] = s2mpjlib('ii','gL4',igt_);
        [it,igt_] = s2mpjlib('ii','gL6',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('(N-2)/3')
            ig = ig_(['OBJ1',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            ig = ig_(['OBJ2',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            ig = ig_(['OBJ3',int2str(I)]);
            pbm.grftype{ig} = 'gL4';
            ig = ig_(['OBJ4',int2str(I)]);
            pbm.grftype{ig} = 'gL6';
        end
        for K=v_('1'):v_('2'):v_('NC')
            v_('K+1') = 1+K;
            ig = ig_(['C',int2str(K)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EA',int2str(K)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['EB',int2str(K)]);
            pbm.grelw{ig}(posel) = 1.;
            ig = ig_(['C',int2str(round(v_('K+1')))]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(round(v_('K+1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               0.0
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COOR2-AY-V-V';
        pbm.objderlvl = 2;
        pb.objderlvl = pbm.objderlvl;
        pbm.conderlvl = [2];
        pb.conderlvl  = pbm.conderlvl;
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb, pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm,pbm.ndigs);
        else
            varargout{1} = pb;
            varargout{2} = pbm;
        end

% **********************
%  SET UP THE FUNCTION *
%  AND RANGE ROUTINES  *
% **********************

    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eC21'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*EV_(1)*EV_(2);
        if(nargout>1)
            g_(1,1) = 2.0*EV_(1)*EV_(2);
            g_(2,1) = EV_(1)*EV_(1);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 2.0*EV_(2);
                H_(1,2) = 2.0*EV_(1);
                H_(2,1) = H_(1,2);
                varargout{3} = H_;
            end
        end

    case 'eS'

        EV_  = varargin{1};
        iel_ = varargin{2};
        U_ = zeros(1,2);
        U_(1,1) = U_(1,1)+1;
        U_(1,2) = U_(1,2)-1;
        IV_(1) = U_(1,:)*EV_;
        SIND = sin(IV_(1));
        varargout{1} = SIND;
        if(nargout>1)
            g_(1,1) = cos(IV_(1));
            varargout{2} = U_.'*g_;
            if(nargout>2)
                H_(1,1) = -SIND;
                varargout{3} = U_.'*H_*U_;
            end
        end

    case 'eC42'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)^4*EV_(2)^2;
        if(nargout>1)
            g_(1,1) = 4.0*EV_(1)^3*EV_(2)^2;
            g_(2,1) = 2.0*EV_(1)^4*EV_(2);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 12.0*EV_(1)^2*EV_(2)^2;
                H_(1,2) = 8.0*EV_(1)^3*EV_(2);
                H_(2,1) = H_(1,2);
                H_(2,2) = 2.0*EV_(1)^4;
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

    case 'gL4'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_^4;
        if(nargout>1)
            g_ = 4.0*GVAR_^3;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 12.0*GVAR_^2;
                varargout{3} = H_;
            end
        end

    case 'gL6'

        GVAR_ = varargin{1};
        igr_  = varargin{2};
        varargout{1} = GVAR_^6;
        if(nargout>1)
            g_ = 6.0*GVAR_^5;
            varargout{2} = g_;
            if(nargout>2)
                H_ = 30.0*GVAR_^4;
                varargout{3} = H_;
            end
        end

    %%%%%%%%%%%%%%%% THE MAIN ACTIONS %%%%%%%%%%%%%%%

    case {'fx','fgx','fgHx','cx','cJx','cJHx','cIx','cIJx','cIJHx','cIJxv','fHxv',...
          'cJxv','cJtxv','cIJtxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy',...
          'LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [0,0];
            [varargout{1:max(1,nargout)}] = s2mpjlib(action,pbm,varargin{:});
        else
            disp(['ERROR: please run ',name,' with action = setup'])
            [varargout{1:nargout}] = deal(NaN);
        end

    otherwise
        disp([' ERROR: action ',action,' unavailable for problem ',name,'.m'])
    end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

