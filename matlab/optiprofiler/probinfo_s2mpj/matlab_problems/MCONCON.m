function varargout = MCONCON(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Another small gas network problem.
% 
%    SIF input: Sybille Schachler, Oxford, August 1992.
%               minor correction by Ph. Shott, Jan 1995.
% 
%    classification = 'C-CLOI2-MN-15-11'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'MCONCON';

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
        v_('N') = 7;
        v_('M') = 4;
        v_('M+1') = 1+v_('M');
        v_('1') = 1;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for I=v_('1'):v_('M')
            [iv,ix_] = s2mpjlib('ii',['P',int2str(I)],ix_);
            pb.xnames{iv} = ['P',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['Q',int2str(I)],ix_);
            pb.xnames{iv} = ['Q',int2str(I)];
            [iv,ix_] = s2mpjlib('ii',['F',int2str(I)],ix_);
            pb.xnames{iv} = ['F',int2str(I)];
        end
        for I=v_('M+1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['P',int2str(I)],ix_);
            pb.xnames{iv} = ['P',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for I=v_('1'):v_('N')
            [ig,ig_] = s2mpjlib('ii','OBJECT',ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['P',int2str(I)]);
            valA(end+1) = -1.0;
        end
        for I=v_('1'):v_('M')
            [ig,ig_] = s2mpjlib('ii',['PAN',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['PAN',int2str(I)];
        end
        [ig,ig_] = s2mpjlib('ii','MBAL1',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'MBAL1';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Q1');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('F3');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','MBAL2',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'MBAL2';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Q1');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('F1');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','MBAL3',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'MBAL3';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Q2');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('F1');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','MBAL4',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'MBAL4';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Q2');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Q3');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','MBAL5',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'MBAL5';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Q3');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('F2');
        valA(end+1) = -1.0;
        [ig,ig_] = s2mpjlib('ii','MBAL6',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'MBAL6';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Q4');
        valA(end+1) = 1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('F2');
        valA(end+1) = 1.0;
        [ig,ig_] = s2mpjlib('ii','MBAL7',ig_);
        gtype{ig}  = '==';
        cnames{ig} = 'MBAL7';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('Q4');
        valA(end+1) = -1.0;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('F4');
        valA(end+1) = -1.0;
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
        v_('DEMAND') = -1000.0;
        pbm.gconst(ig_('MBAL4')) = v_('DEMAND');
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        v_('PMAX1') = 914.73;
        v_('PMAX2') = 904.73;
        pb.xupper(ix_('P3')) = v_('PMAX2');
        pb.xupper(ix_('P5')) = v_('PMAX2');
        pb.xupper(ix_('P1')) = v_('PMAX1');
        pb.xupper(ix_('P7')) = v_('PMAX1');
        pb.xupper(ix_('F4')) = 400.0;
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        for I=v_('1'):v_('N')
            pb.x0(ix_(['P',int2str(I)]),1) = 965.0;
        end
        pb.x0(ix_('Q1'),1) = 100.0;
        pb.x0(ix_('Q2'),1) = 100.0;
        pb.x0(ix_('Q3'),1) = -100.0;
        pb.x0(ix_('Q4'),1) = -100.0;
        pb.x0(ix_('F1'),1) = 1000.0;
        pb.x0(ix_('F2'),1) = 1000.0;
        pb.x0(ix_('F3'),1) = 1000.0;
        pb.x0(ix_('F4'),1) = 1000.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eSQR',iet_);
        elftv{it}{1} = 'X';
        [it,iet_] = s2mpjlib( 'ii', 'eFORQ',iet_);
        elftv{it}{1} = 'Y';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('N')
            ename = ['PSQ',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eSQR';
            ielftype(ie) = iet_('eSQR');
            vname = ['P',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        for I=v_('1'):v_('M')
            ename = ['QTO',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eFORQ';
            ielftype(ie) = iet_('eFORQ');
            vname = ['Q',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        v_('K') = -0.597053452;
        ig = ig_('PAN1');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('PSQ1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('PSQ2');
        pbm.grelw{ig}(posel) = -1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('QTO1');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('K');
        ig = ig_('PAN2');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('PSQ3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('PSQ4');
        pbm.grelw{ig}(posel) = -1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('QTO2');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('K');
        ig = ig_('PAN3');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('PSQ4');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('PSQ5');
        pbm.grelw{ig}(posel) = -1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('QTO3');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('K');
        ig = ig_('PAN4');
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('PSQ6');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = 1.0;
        posel = posel+1;
        pbm.grelt{ig}(posel) = ie_('PSQ7');
        pbm.grelw{ig}(posel) = -1.0;
        posel = length(pbm.grelt{ig})+1;
        pbm.grelt{ig}(posel) = ie_('QTO4');
        nlc = union(nlc,ig);
        pbm.grelw{ig}(posel) = v_('K');
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CLOI2-MN-15-11';
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


    %%%%%%%%%%%%%%%% NONLINEAR ELEMENTS %%%%%%%%%%%%%%%

    case 'eSQR'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*abs(EV_(1));
        PO = EV_(1)>0.0;
        if(PO)
            GO = 2*EV_(1);
        end
        if(~PO)
            GO = -2*EV_(1);
        end
        if(PO)
            HO = 2;
        end
        if(~PO)
            HO = -2;
        end
        if(nargout>1)
            g_(1,1) = GO;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = HO;
                varargout{3} = H_;
            end
        end

    case 'eFORQ'

        EV_  = varargin{1};
        iel_ = varargin{2};
        varargout{1} = EV_(1)*abs(EV_(1))^0.8539;
        POS = EV_(1)>0.0;
        if(POS)
            GG = 1.8539*EV_(1)^0.8539;
        end
        if(~POS)
            GG = 1.8539*abs(EV_(1))^0.8539;
        end
        if(POS)
            HH = 1.8539*0.8539*EV_(1)^(-0.1461);
        end
        if(~POS)
            HH = -1.8539*0.8539*abs(EV_(1))^(-0.1461);
        end
        if(nargout>1)
            g_(1,1) = GG;
            varargout{2} = g_;
            if(nargout>2)
                H_(1,1) = HH;
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

