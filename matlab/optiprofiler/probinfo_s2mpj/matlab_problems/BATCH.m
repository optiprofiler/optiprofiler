function varargout = BATCH(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : BATCH
%    *********
% 
%    Source: Optimal Design of Multiproduct Batch Plant
%    G.R. Kocis & I.E. Grossmann,
%    "Global OPtimization of Nonconvex Mixed Integer Nonlinear Programmming
%     (MINLP) problems in Process Synthesis", Indust. Engng. Chem. Res.,
%    No. 27, pp 1407--1421, 1988.
% 
%    SIF input: S. Leyffer, October 1997
% 
%    classification = 'C-COOR2-AN-46-73'
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 17 X 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'BATCH';

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
        v_  = containers.Map('KeyType','char', 'ValueType', 'double');
        ix_ = containers.Map('KeyType','char', 'ValueType', 'double');
        ig_ = containers.Map('KeyType','char', 'ValueType', 'double');
        v_('1') = 1;
        v_('M') = 6;
        v_('N') = 5;
        v_('NU') = 4;
        v_('LOGNU') = log(4.0);
        v_('VL') = log(300.0);
        v_('VU') = log(3000.0);
        v_('H') = 6000.0;
        v_('TLO1') = 0.729961;
        v_('TLO2') = 0.530628;
        v_('TLO3') = 1.09024;
        v_('TLO4') = -0.133531;
        v_('TLO5') = 0.0487901;
        v_('TUP1') = 2.11626;
        v_('TUP2') = 1.91626;
        v_('TUP3') = 2.47654;
        v_('TUP4') = 1.25276;
        v_('TUP5') = 1.43508;
        v_('BLO1') = 4.45966;
        v_('BLO2') = 3.74950;
        v_('BLO3') = 4.49144;
        v_('BLO4') = 3.14988;
        v_('BLO5') = 3.04452;
        v_('BUP1') = 397.747;
        v_('BUP2') = 882.353;
        v_('BUP3') = 833.333;
        v_('BUP4') = 638.298;
        v_('BUP5') = 666.667;
        v_('Q1') = 250000.0;
        v_('Q2') = 150000.0;
        v_('Q3') = 180000.0;
        v_('Q4') = 160000.0;
        v_('Q5') = 120000.0;
        v_('LOGI1') = log(1.0);
        v_('LOGI2') = log(2.0);
        v_('LOGI3') = log(3.0);
        v_('LOGI4') = log(4.0);
        v_('S1,1') = log(7.9);
        v_('S2,1') = log(0.7);
        v_('S3,1') = log(0.7);
        v_('S4,1') = log(4.7);
        v_('S5,1') = log(1.2);
        v_('S1,2') = log(2.0);
        v_('S2,2') = log(0.8);
        v_('S3,2') = log(2.6);
        v_('S4,2') = log(2.3);
        v_('S5,2') = log(3.6);
        v_('S1,3') = log(5.2);
        v_('S2,3') = log(0.9);
        v_('S3,3') = log(1.6);
        v_('S4,3') = log(1.6);
        v_('S5,3') = log(2.4);
        v_('S1,4') = log(4.9);
        v_('S2,4') = log(3.4);
        v_('S3,4') = log(3.6);
        v_('S4,4') = log(2.7);
        v_('S5,4') = log(4.5);
        v_('S1,5') = log(6.1);
        v_('S2,5') = log(2.1);
        v_('S3,5') = log(3.2);
        v_('S4,5') = log(1.2);
        v_('S5,5') = log(1.6);
        v_('S1,6') = log(4.2);
        v_('S2,6') = log(2.5);
        v_('S3,6') = log(2.9);
        v_('S4,6') = log(2.5);
        v_('S5,6') = log(2.1);
        v_('T1,1') = log(6.4);
        v_('T2,1') = log(6.8);
        v_('T3,1') = log(1.0);
        v_('T4,1') = log(3.2);
        v_('T5,1') = log(2.1);
        v_('T1,2') = log(4.7);
        v_('T2,2') = log(6.4);
        v_('T3,2') = log(6.3);
        v_('T4,2') = log(3.0);
        v_('T5,2') = log(2.5);
        v_('T1,3') = log(8.3);
        v_('T2,3') = log(6.5);
        v_('T3,3') = log(5.4);
        v_('T4,3') = log(3.5);
        v_('T5,3') = log(4.2);
        v_('T1,4') = log(3.9);
        v_('T2,4') = log(4.4);
        v_('T3,4') = log(11.9);
        v_('T4,4') = log(3.3);
        v_('T5,4') = log(3.6);
        v_('T1,5') = log(2.1);
        v_('T2,5') = log(2.3);
        v_('T3,5') = log(5.7);
        v_('T4,5') = log(2.8);
        v_('T5,5') = log(3.7);
        v_('T1,6') = log(1.2);
        v_('T2,6') = log(3.2);
        v_('T3,6') = log(6.2);
        v_('T4,6') = log(3.4);
        v_('T5,6') = log(2.2);
        for J=v_('1'):v_('M')
            v_(['ALPHA',int2str(J)]) = 250.0;
            v_(['BETA',int2str(J)]) = 0.6;
        end
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for J=v_('1'):v_('M')
            [iv,ix_] = s2mpjlib('ii',['N',int2str(J)],ix_);
            pb.xnames{iv} = ['N',int2str(J)];
        end
        for J=v_('1'):v_('M')
            [iv,ix_] = s2mpjlib('ii',['V',int2str(J)],ix_);
            pb.xnames{iv} = ['V',int2str(J)];
        end
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['B',int2str(I)],ix_);
            pb.xnames{iv} = ['B',int2str(I)];
        end
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['TL',int2str(I)],ix_);
            pb.xnames{iv} = ['TL',int2str(I)];
        end
        for J=v_('1'):v_('M')
            for K=v_('1'):v_('NU')
                [iv,ix_] = s2mpjlib('ii',['Y',int2str(K),',',int2str(J)],ix_);
                pb.xnames{iv} = ['Y',int2str(K),',',int2str(J)];
            end
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        [ig,ig_] = s2mpjlib('ii','COST',ig_);
        gtype{ig} = '<>';
        for I=v_('1'):v_('N')
            for J=v_('1'):v_('M')
                [ig,ig_] = s2mpjlib('ii',['VOL',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '>=';
                cnames{ig} = ['VOL',int2str(I),',',int2str(J)];
                iv = ix_(['V',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                iv = ix_(['B',int2str(I)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = -1.0;
                end
            end
        end
        for I=v_('1'):v_('N')
            for J=v_('1'):v_('M')
                [ig,ig_] = s2mpjlib('ii',['CYCL',int2str(I),',',int2str(J)],ig_);
                gtype{ig}  = '>=';
                cnames{ig} = ['CYCL',int2str(I),',',int2str(J)];
                iv = ix_(['N',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
                iv = ix_(['TL',int2str(I)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
            end
        end
        [ig,ig_] = s2mpjlib('ii','HORIZON',ig_);
        gtype{ig}  = '<=';
        cnames{ig} = 'HORIZON';
        for J=v_('1'):v_('M')
            for K=v_('1'):v_('NU')
                [ig,ig_] = s2mpjlib('ii',['NPAR',int2str(J)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['NPAR',int2str(J)];
                iv = ix_(['Y',int2str(K),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = v_(['LOGI',int2str(K)])+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = v_(['LOGI',int2str(K)]);
                end
            end
            [ig,ig_] = s2mpjlib('ii',['NPAR',int2str(J)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['NPAR',int2str(J)];
            iv = ix_(['N',int2str(J)]);
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = -1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = -1.0;
            end
        end
        for J=v_('1'):v_('M')
            for K=v_('1'):v_('NU')
                [ig,ig_] = s2mpjlib('ii',['SOS1',int2str(J)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['SOS1',int2str(J)];
                iv = ix_(['Y',int2str(K),',',int2str(J)]);
                if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                    pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
                else
                    pbm.A(ig,iv) = 1.0;
                end
            end
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
        for I=v_('1'):v_('N')
            for J=v_('1'):v_('M')
                pbm.gconst(ig_(['VOL',int2str(I),',',int2str(J)])) =...
                      v_(['S',int2str(I),',',int2str(J)]);
            end
        end
        for I=v_('1'):v_('N')
            for J=v_('1'):v_('M')
                pbm.gconst(ig_(['CYCL',int2str(I),',',int2str(J)])) =...
                      v_(['T',int2str(I),',',int2str(J)]);
            end
        end
        pbm.gconst(ig_('HORIZON')) = v_('H');
        for J=v_('1'):v_('M')
            pbm.gconst(ig_(['SOS1',int2str(J)])) = 1.0;
        end
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        for J=v_('1'):v_('M')
            pb.xupper(ix_(['N',int2str(J)])) = v_('LOGNU');
            pb.xlower(ix_(['V',int2str(J)]),1) = v_('VL');
            pb.xupper(ix_(['V',int2str(J)])) = v_('VU');
        end
        for I=v_('1'):v_('N')
            pb.xlower(ix_(['B',int2str(I)]),1) = v_(['BLO',int2str(I)]);
            pb.xupper(ix_(['B',int2str(I)])) = v_(['BUP',int2str(I)]);
            pb.xlower(ix_(['TL',int2str(I)]),1) = v_(['TLO',int2str(I)]);
            pb.xupper(ix_(['TL',int2str(I)])) = v_(['TUP',int2str(I)]);
        end
        for J=v_('1'):v_('M')
            for K=v_('1'):v_('NU')
                pb.xupper(ix_(['Y',int2str(K),',',int2str(J)])) = 1.0;
            end
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eEXPXAY',iet_);
        elftv{it}{1} = 'X';
        elftv{it}{2} = 'Y';
        elftp{it}{1} = 'A';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for J=v_('1'):v_('M')
            ename = ['EXPO',int2str(J)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eEXPXAY';
            ielftype(ie) = iet_('eEXPXAY');
            vname = ['N',int2str(J)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['V',int2str(J)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('A',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['BETA',int2str(J)]);
        end
        for I=v_('1'):v_('M')
            ename = ['EXPC',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eEXPXAY';
            ielftype(ie) = iet_('eEXPXAY');
            vname = ['TL',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['B',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('Y',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('A',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = -1.0;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for J=v_('1'):v_('M')
            ig = ig_('COST');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EXPO',int2str(J)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_(['ALPHA',int2str(J)]);
        end
        for I=v_('1'):v_('N')
            ig = ig_('HORIZON');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EXPC',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_(['Q',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Solution
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(1:pb.nle) = -Inf*ones(pb.nle,1);
        pb.cupper(1:pb.nle) = zeros(pb.nle,1);
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.clower(pb.nle+pb.neq+1:pb.m) = zeros(pb.nge,1);
        pb.cupper(1:pb.nge) = +Inf*ones(pb.nge,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-COOR2-AN-46-73';
        pb.x0          = zeros(pb.n,1);
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

    case 'eEXPXAY'

        EV_  = varargin{1};
        iel_ = varargin{2};
        FVALUE = exp(EV_(1)+pbm.elpar{iel_}(1)*EV_(2));
        GYVALU = pbm.elpar{iel_}(1)*FVALUE;
        varargout{1} = FVALUE;
        if(nargout>1)
            g_(1,1) = FVALUE;
            g_(2,1) = GYVALU;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = FVALUE;
                H_(1,2) = GYVALU;
                H_(2,1) = H_(1,2);
                H_(2,2) = pbm.elpar{iel_}(1)*GYVALU;
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

