function varargout = ELATTAR(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : ELATTAR
%    *********
% 
%    A nonlinear minmax problem in six variables.
% 
%    The problem is nonconvex and has several local minima.
% 
%    Source: 
%    R.A. El-Attar, M. Vidyasagar and S.R.K. Dutta,
%    "An algorithm for l_1-approximation",
%    SINUM 16, pp.70-86, 1979.
% 
%    SIF input: Ph. Toint, Nov 1993.
% 
%    classification = 'C-CLOR2-AN-7-102'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'ELATTAR';

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
        v_('1') = 1;
        v_('6') = 6;
        v_('51') = 51;
        v_('T') = 0.0;
        for I=v_('1'):v_('51')
            v_(['T',int2str(I)]) = v_('T');
            v_('T') = 0.1+v_('T');
            v_('ETI') = exp(v_(['T',int2str(I)]));
            v_(['Y',int2str(I)]) = 0.5*v_('ETI');
            v_('-2TI') = -2.0*v_(['T',int2str(I)]);
            v_('E-2TI') = exp(v_('-2TI'));
            v_(['Y',int2str(I)]) = v_(['Y',int2str(I)])-v_('E-2TI');
            v_('-3TI') = -3.0*v_(['T',int2str(I)]);
            v_('E-3TI') = exp(v_('-3TI'));
            v_('E-3TI/2') = 0.5*v_('E-3TI');
            v_(['Y',int2str(I)]) = v_(['Y',int2str(I)])+v_('E-3TI/2');
            v_('-3TI/2') = 0.5*v_('-3TI');
            v_('E-3TI/2') = exp(v_('-3TI/2'));
            v_('7TI') = 7.0*v_(['T',int2str(I)]);
            v_('S7TI') = sin(v_('7TI'));
            v_('TT') = v_('E-3TI/2')*v_('S7TI');
            v_('TT') = 1.5*v_('TT');
            v_(['Y',int2str(I)]) = v_(['Y',int2str(I)])+v_('TT');
            v_('5TI') = 5.0*v_(['T',int2str(I)]);
            v_('-5TI/2') = -0.5*v_('5TI');
            v_('E-5TI/2') = exp(v_('-5TI/2'));
            v_('S5TI') = sin(v_('5TI'));
            v_('TT') = v_('E-5TI/2')*v_('S5TI');
            v_(['Y',int2str(I)]) = v_(['Y',int2str(I)])+v_('TT');
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('6')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(I)],ix_);
            pb.xnames{iv} = ['X',int2str(I)];
        end
        [iv,ix_] = s2mpjlib('ii','U',ix_);
        pb.xnames{iv} = 'U';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('U');
        valA(end+1) = 1.0;
        for I=v_('1'):v_('51')
            [ig,ig_] = s2mpjlib('ii',['F',int2str(I)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['F',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_('U');
            valA(end+1) = -1.0;
            [ig,ig_] = s2mpjlib('ii',['MF',int2str(I)],ig_);
            gtype{ig}  = '<=';
            cnames{ig} = ['MF',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_('U');
            valA(end+1) = -1.0;
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
        for I=v_('1'):v_('51')
            pbm.gconst(ig_(['F',int2str(I)])) = v_(['Y',int2str(I)]);
            v_(['-Y',int2str(I)]) = -1.0*v_(['Y',int2str(I)]);
            pbm.gconst(ig_(['MF',int2str(I)])) = v_(['-Y',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'X1'))
            pb.x0(ix_('X1'),1) = -2.0;
        else
            pb.y0(find(pbm.congrps==ig_('X1')),1) = -2.0;
        end
        if(isKey(ix_,'X2'))
            pb.x0(ix_('X2'),1) = -2.0;
        else
            pb.y0(find(pbm.congrps==ig_('X2')),1) = -2.0;
        end
        if(isKey(ix_,'X3'))
            pb.x0(ix_('X3'),1) = 7.0;
        else
            pb.y0(find(pbm.congrps==ig_('X3')),1) = 7.0;
        end
        if(isKey(ix_,'X5'))
            pb.x0(ix_('X5'),1) = -2.0;
        else
            pb.y0(find(pbm.congrps==ig_('X5')),1) = -2.0;
        end
        if(isKey(ix_,'X6'))
            pb.x0(ix_('X6'),1) = 1.0;
        else
            pb.y0(find(pbm.congrps==ig_('X6')),1) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eET1',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftv{it}{3} = 'V3';
        elftv{it}{4} = 'V4';
        elftp{it}{1} = 'T';
        [it,iet_] = s2mpjlib( 'ii', 'eET2',iet_);
        elftv{it}{1} = 'V5';
        elftv{it}{2} = 'V6';
        elftp{it}{1} = 'T';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('51')
            ename = ['EL1',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eET1';
            ielftype(ie) = iet_('eET1');
            vname = 'X1';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X4';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V4',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('T',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['T',int2str(I)]);
            ename = ['EL2',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eET2';
            ielftype(ie) = iet_('eET2');
            vname = 'X5';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V5',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X6';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V6',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('T',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['T',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('51')
            ig = ig_(['F',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EL1',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.0;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['EL2',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.0;
            ig = ig_(['MF',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EL1',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = -1.0;
            posel = posel+1;
            pbm.grelt{ig}(posel) = ie_(['EL2',int2str(I)]);
            pbm.grelw{ig}(posel) = -1.0;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution       
% LO SOLTN               0.1427066255
% LO SOLTN               74.206179244
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CLOR2-AN-7-102';
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

    case 'eET1'

        EV_  = varargin{1};
        iel_ = varargin{2};
        A = -EV_(2)*pbm.elpar{iel_}(1);
        B = EV_(3)*pbm.elpar{iel_}(1)+EV_(4);
        EA = exp(A);
        CB = cos(B);
        SB = sin(B);
        EACB = EA*CB;
        EASB = EA*SB;
        V1EACB = EV_(1)*EACB;
        V1EASB = EV_(1)*EASB;
        varargout{1} = V1EACB;
        if(nargout>1)
            g_(1,1) = EACB;
            g_(2,1) = -pbm.elpar{iel_}(1)*V1EACB;
            g_(3,1) = -pbm.elpar{iel_}(1)*V1EASB;
            g_(4,1) = -V1EASB;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(4,4);
                H_(1,2) = -pbm.elpar{iel_}(1)*EACB;
                H_(2,1) = H_(1,2);
                H_(1,3) = -pbm.elpar{iel_}(1)*EASB;
                H_(3,1) = H_(1,3);
                H_(1,4) = -EASB;
                H_(4,1) = H_(1,4);
                H_(2,2) = pbm.elpar{iel_}(1)*pbm.elpar{iel_}(1)*V1EACB;
                H_(2,3) = pbm.elpar{iel_}(1)*pbm.elpar{iel_}(1)*V1EASB;
                H_(3,2) = H_(2,3);
                H_(2,4) = pbm.elpar{iel_}(1)*V1EASB;
                H_(4,2) = H_(2,4);
                H_(3,3) = -pbm.elpar{iel_}(1)*pbm.elpar{iel_}(1)*V1EACB;
                H_(3,4) = -pbm.elpar{iel_}(1)*V1EACB;
                H_(4,3) = H_(3,4);
                H_(4,4) = -V1EACB;
                varargout{3} = H_;
            end
        end

    case 'eET2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        A = -EV_(2)*pbm.elpar{iel_}(1);
        EA = exp(A);
        B = EV_(1)*EA;
        varargout{1} = B;
        if(nargout>1)
            g_(1,1) = EA;
            g_(2,1) = -pbm.elpar{iel_}(1)*B;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,2) = -pbm.elpar{iel_}(1)*EA;
                H_(2,1) = H_(1,2);
                H_(2,2) = pbm.elpar{iel_}(1)*pbm.elpar{iel_}(1)*B;
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

