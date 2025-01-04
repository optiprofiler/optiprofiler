function varargout = HIMMELBK(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : HIMMELBK
%    *********
% 
%    A blending problem for multi-component mixtures, by Paviani.
%    It has a linear objective and linear and nonlinear constraints.
% 
%    Compared to the problem specified in Himmelblau, the inequality
%    constraints have been removed, because, as stated in this source,
%    they impose that
%    X(1)=X(2)=X(3)=X(7)=X(9)=X(9)=X(13)=X(14)=X(15)=X(19)=X(20)=X(21)=0
%    which is clearly contradictory with the given solution(s).  As
%    there does not seem to be a natural way to correct this statement
%    without knowing more about the original problem, the troublesome
%    constraints have been removed.
% 
%    Source: from problem 20 in
%    D.H. Himmelblau,
%    "Applied nonlinear programming",
%    McGraw-Hill, New-York, 1972.
% 
%    SIF input: Ph. Toint, March 1991.
% 
%    classification = 'C-CLOR2-MN-24-14'
% 
%    Problem data
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'HIMMELBK';

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
        v_('F') = 142.22471;
        v_('B1') = 44.094;
        v_('B2') = 58.12;
        v_('B3') = 58.12;
        v_('B4') = 137.4;
        v_('B5') = 120.9;
        v_('B6') = 170.9;
        v_('B7') = 62.501;
        v_('B8') = 84.94;
        v_('B9') = 133.425;
        v_('B10') = 82.507;
        v_('B11') = 46.07;
        v_('B12') = 60.097;
        v_('B13') = 44.094;
        v_('B14') = 58.12;
        v_('B15') = 58.12;
        v_('B16') = 137.4;
        v_('B17') = 120.9;
        v_('B18') = 170.9;
        v_('B19') = 62.501;
        v_('B20') = 84.94;
        v_('B21') = 133.425;
        v_('B22') = 82.507;
        v_('B23') = 46.07;
        v_('B24') = 60.097;
        v_('C1') = 123.7;
        v_('C2') = 31.7;
        v_('C3') = 45.7;
        v_('C4') = 14.7;
        v_('C5') = 84.7;
        v_('C6') = 27.7;
        v_('C7') = 49.7;
        v_('C8') = 7.1;
        v_('C9') = 2.1;
        v_('C10') = 17.7;
        v_('C11') = 0.85;
        v_('C12') = 0.64;
        v_('D1') = 123.7;
        v_('D2') = 31.7;
        v_('D3') = 45.7;
        v_('D4') = 14.7;
        v_('D5') = 84.7;
        v_('D6') = 27.7;
        v_('D7') = 49.7;
        v_('D8') = 7.1;
        v_('D9') = 2.1;
        v_('D10') = 17.7;
        v_('D11') = 0.85;
        v_('D12') = 0.64;
        v_('1') = 1;
        v_('12') = 12;
        v_('13') = 13;
        v_('24') = 24;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        for K=v_('1'):v_('24')
            [iv,ix_] = s2mpjlib('ii',['X',int2str(K)],ix_);
            pb.xnames{iv} = ['X',int2str(K)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X1');
        valA(end+1) = 0.0693;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X2');
        valA(end+1) = 0.0577;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X3');
        valA(end+1) = 0.05;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X4');
        valA(end+1) = 0.2;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X5');
        valA(end+1) = 0.26;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X6');
        valA(end+1) = 0.55;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X7');
        valA(end+1) = 0.06;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X8');
        valA(end+1) = 0.1;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X9');
        valA(end+1) = 0.12;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X10');
        valA(end+1) = 0.18;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X11');
        valA(end+1) = 0.1;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X12');
        valA(end+1) = 0.09;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X13');
        valA(end+1) = 0.0693;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X14');
        valA(end+1) = 0.0577;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X15');
        valA(end+1) = 0.05;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X16');
        valA(end+1) = 0.2;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X17');
        valA(end+1) = 0.26;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X18');
        valA(end+1) = 0.55;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X19');
        valA(end+1) = 0.06;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X20');
        valA(end+1) = 0.1;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X21');
        valA(end+1) = 0.12;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X22');
        valA(end+1) = 0.18;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X23');
        valA(end+1) = 0.1;
        irA(end+1)  = ig;
        icA(end+1)  = ix_('X24');
        valA(end+1) = 0.09;
        for I=v_('1'):v_('12')
            [ig,ig_] = s2mpjlib('ii',['CA',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['CA',int2str(I)];
        end
        for I=v_('1'):v_('24')
            [ig,ig_] = s2mpjlib('ii','CA13',ig_);
            gtype{ig}  = '==';
            cnames{ig} = 'CA13';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I)]);
            valA(end+1) = 1.0;
        end
        for I=v_('1'):v_('12')
            v_('I+12') = 12+I;
            v_('1/DI') = 1.0/v_(['D',int2str(I)]);
            [ig,ig_] = s2mpjlib('ii','CA14',ig_);
            gtype{ig}  = '==';
            cnames{ig} = 'CA14';
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(I)]);
            valA(end+1) = v_('1/DI');
            v_('F/BI+12') = v_('F')/v_(['B',int2str(round(v_('I+12')))]);
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['X',int2str(round(v_('I+12')))]);
            valA(end+1) = v_('F/BI+12');
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
        pbm.gconst(ig_('CA13')) = 1.0;
        pbm.gconst(ig_('CA14')) = 1.671;
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.04*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'en2PR',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        for I=v_('1'):v_('12')
            v_('I+12') = 12+I;
            for J=v_('1'):v_('12')
                ename = ['E',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'en2PR';
                ielftype(ie) = iet_('en2PR');
                vname = ['X',int2str(round(v_('I+12')))];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.04);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.04);
                posev = find(strcmp('Y',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
            for J=v_('13'):v_('24')
                ename = ['E',int2str(I),',',int2str(J)];
                [ie,ie_] = s2mpjlib('ii',ename,ie_);
                pbm.elftype{ie} = 'en2PR';
                ielftype(ie) = iet_('en2PR');
                vname = ['X',int2str(I)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.04);
                posev = find(strcmp('X',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
                vname = ['X',int2str(J)];
                [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.04);
                posev = find(strcmp('Y',elftv{ielftype(ie)}));
                pbm.elvar{ie}(posev) = iv;
            end
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('12')
            v_('I+12') = 12+I;
            for J=v_('1'):v_('12')
                v_('BI/BJ') = v_(['B',int2str(I)])/v_(['B',int2str(J)]);
                v_('40BI/BJ') = 40.0*v_('BI/BJ');
                ig = ig_(['CA',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('40BI/BJ');
            end
            for J=v_('13'):v_('24')
                v_('B+/BJ') = v_(['B',int2str(round(v_('I+12')))])/v_(['B',int2str(J)]);
                v_('CB+/BJ') = v_(['C',int2str(I)])*v_('B+/BJ');
                v_('-CB+/BJ') = -1.0*v_('CB+/BJ');
                ig = ig_(['CA',int2str(I)]);
                posel = length(pbm.grelt{ig})+1;
                pbm.grelt{ig}(posel) = ie_(['E',int2str(I),',',int2str(J)]);
                nlc = union(nlc,ig);
                pbm.grelw{ig}(posel) = v_('-CB+/BJ');
            end
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
% LO SOLTN                0.0893344
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CLOR2-MN-24-14';
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

