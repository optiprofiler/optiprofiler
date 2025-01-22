function varargout = EXPFITB(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : EXPFITB
%    *********
% 
%    One sided rational approximation to the exponential function, as
%    described by Powell.
% 
%    Source:
%    M.J.D. Powell,
%    "A tolerant algorithm for linearly constrained optimization
%    calculations"'
%    Mathematical Programming 45(3), pp.561--562, 1989.
% 
%    SDIF input: Ph. Toint and N. Gould, May 1990.
% 
%    classification = 'C-COLR2-AN-5-102'
% 
%    Number of fitting points
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'EXPFITB';

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
        v_('R') = 51;
        v_('1') = 1;
        v_('5.0') = 5.0;
        v_('R-1') = -1+v_('R');
        v_('RR-1') = v_('R-1');
        v_('5/R-1') = v_('5.0')/v_('RR-1');
        for I=v_('1'):v_('R')
            v_('I-1') = -1+I;
            v_('RI-1') = v_('I-1');
            v_(['T',int2str(I)]) = v_('RI-1')*v_('5/R-1');
            v_(['ET',int2str(I)]) = exp(v_(['T',int2str(I)]));
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        [iv,ix_] = s2mpjlib('ii','P0',ix_);
        pb.xnames{iv} = 'P0';
        [iv,ix_] = s2mpjlib('ii','P1',ix_);
        pb.xnames{iv} = 'P1';
        [iv,ix_] = s2mpjlib('ii','P2',ix_);
        pb.xnames{iv} = 'P2';
        [iv,ix_] = s2mpjlib('ii','Q1',ix_);
        pb.xnames{iv} = 'Q1';
        [iv,ix_] = s2mpjlib('ii','Q2',ix_);
        pb.xnames{iv} = 'Q2';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        [ig,ig_] = s2mpjlib('ii','OBJ',ig_);
        gtype{ig} = '<>';
        for I=v_('1'):v_('R')
            v_('TM5') = -5.0+v_(['T',int2str(I)]);
            v_('TM5SQ') = v_('TM5')*v_('TM5');
            v_('QC1') = v_('TM5')*v_(['ET',int2str(I)]);
            v_('QC2') = v_('TM5SQ')*v_(['ET',int2str(I)]);
            v_('-QC1') = -1.0*v_('QC1');
            v_('-QC2') = -1.0*v_('QC2');
            v_('2T') = v_(['T',int2str(I)])*v_(['T',int2str(I)]);
            [ig,ig_] = s2mpjlib('ii',['C',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['C',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_('P0');
            valA(end+1) = 1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_('P1');
            valA(end+1) = v_(['T',int2str(I)]);
            irA(end+1)  = ig;
            icA(end+1)  = ix_('P2');
            valA(end+1) = v_('2T');
            irA(end+1)  = ig;
            icA(end+1)  = ix_('Q1');
            valA(end+1) = v_('-QC1');
            irA(end+1)  = ig;
            icA(end+1)  = ix_('Q2');
            valA(end+1) = v_('-QC2');
            [ig,ig_] = s2mpjlib('ii',['B',int2str(I)],ig_);
            gtype{ig}  = '>=';
            cnames{ig} = ['B',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_('Q1');
            valA(end+1) = v_('TM5');
            irA(end+1)  = ig;
            icA(end+1)  = ix_('Q2');
            valA(end+1) = v_('TM5SQ');
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
        for I=v_('1'):v_('R')
            pbm.gconst(ig_(['C',int2str(I)])) = v_(['ET',int2str(I)]);
            pbm.gconst(ig_(['B',int2str(I)])) = -0.99999;
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'P0'))
            pb.x0(ix_('P0'),1) = 1.0;
        else
            pb.y0(find(pbm.congrps==ig_('P0')),1) = 1.0;
        end
        if(isKey(ix_,'P1'))
            pb.x0(ix_('P1'),1) = 1.0;
        else
            pb.y0(find(pbm.congrps==ig_('P1')),1) = 1.0;
        end
        if(isKey(ix_,'P2'))
            pb.x0(ix_('P2'),1) = 6.0;
        else
            pb.y0(find(pbm.congrps==ig_('P2')),1) = 6.0;
        end
        if(isKey(ix_,'Q1'))
            pb.x0(ix_('Q1'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('Q1')),1) = 0.0;
        end
        if(isKey(ix_,'Q2'))
            pb.x0(ix_('Q2'),1) = 0.0;
        else
            pb.y0(find(pbm.congrps==ig_('Q2')),1) = 0.0;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eFIT',iet_);
        elftv{it}{1} = 'P0';
        elftv{it}{2} = 'P1';
        elftv{it}{3} = 'P2';
        elftv{it}{4} = 'Q1';
        elftv{it}{5} = 'Q2';
        elftp{it}{1} = 'T';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('R')
            ename = ['F',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eFIT';
            ielftype(ie) = iet_('eFIT');
            vname = 'P0';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('P0',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'P1';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('P1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'P2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('P2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'Q1';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Q1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'Q2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Q2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('T',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['T',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('R')
            ig = ig_('OBJ');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['F',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COLR2-AN-5-102';
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

    case 'eFIT'

        EV_  = varargin{1};
        iel_ = varargin{2};
        TM5 = pbm.elpar{iel_}(1)-5.0;
        TM5SQ = TM5*TM5;
        T2 = pbm.elpar{iel_}(1)*pbm.elpar{iel_}(1);
        ET = exp(pbm.elpar{iel_}(1));
        QT = 1.0+EV_(4)*TM5+EV_(5)*TM5SQ;
        ETQT = ET*QT;
        ETQT2 = ETQT*QT;
        ETQT3 = ETQT2*QT;
        PT = EV_(1)+EV_(2)*pbm.elpar{iel_}(1)+EV_(3)*T2;
        F = PT/ETQT-1.0;
        TWOF = F+F;
        DFDP0 = 1.0/ETQT;
        DFDP1 = pbm.elpar{iel_}(1)/ETQT;
        DFDP2 = T2/ETQT;
        DFDQ1 = -PT*TM5/ETQT2;
        DFDQ2 = -PT*TM5SQ/ETQT2;
        D2P0Q1 = -TM5/ETQT2;
        D2P0Q2 = -TM5SQ/ETQT2;
        D2P1Q1 = -pbm.elpar{iel_}(1)*TM5/ETQT2;
        D2P1Q2 = -pbm.elpar{iel_}(1)*TM5SQ/ETQT2;
        D2P2Q1 = -T2*TM5/ETQT2;
        D2P2Q2 = -T2*TM5SQ/ETQT2;
        D2Q1Q1 = 2.0*PT*TM5SQ/ETQT3;
        D2Q1Q2 = 2.0*PT*TM5SQ*TM5/ETQT3;
        D2Q2Q2 = 2.0*PT*TM5SQ*TM5SQ/ETQT3;
        varargout{1} = F*F;
        if(nargout>1)
            g_(1,1) = TWOF*DFDP0;
            g_(2,1) = TWOF*DFDP1;
            g_(3,1) = TWOF*DFDP2;
            g_(4,1) = TWOF*DFDQ1;
            g_(5,1) = TWOF*DFDQ2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(5,5);
                H_(1,1) = 2.0*DFDP0*DFDP0;
                H_(1,2) = 2.0*DFDP0*DFDP1;
                H_(2,1) = H_(1,2);
                H_(1,3) = 2.0*DFDP0*DFDP2;
                H_(3,1) = H_(1,3);
                H_(2,2) = 2.0*DFDP1*DFDP1;
                H_(2,3) = 2.0*DFDP1*DFDP2;
                H_(3,2) = H_(2,3);
                H_(3,3) = 2.0*DFDP2*DFDP2;
                H_(1,4) = TWOF*D2P0Q1+2.0*DFDP0*DFDQ1;
                H_(4,1) = H_(1,4);
                H_(1,5) = TWOF*D2P0Q2+2.0*DFDP0*DFDQ2;
                H_(5,1) = H_(1,5);
                H_(2,4) = TWOF*D2P1Q1+2.0*DFDP1*DFDQ1;
                H_(4,2) = H_(2,4);
                H_(2,5) = TWOF*D2P1Q2+2.0*DFDP1*DFDQ2;
                H_(5,2) = H_(2,5);
                H_(3,4) = TWOF*D2P2Q1+2.0*DFDP2*DFDQ1;
                H_(4,3) = H_(3,4);
                H_(3,5) = TWOF*D2P2Q2+2.0*DFDP2*DFDQ2;
                H_(5,3) = H_(3,5);
                H_(4,4) = TWOF*D2Q1Q1+2.0*DFDQ1*DFDQ1;
                H_(4,5) = TWOF*D2Q1Q2+2.0*DFDQ1*DFDQ2;
                H_(5,4) = H_(4,5);
                H_(5,5) = TWOF*D2Q2Q2+2.0*DFDQ2*DFDQ2;
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

