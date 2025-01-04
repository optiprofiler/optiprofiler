function varargout = LUBRIFC(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : LUBRIFC
%    *********
% 
%    Corrected version of LUBRIF which contained an error
%    in the definition of the Reynold's equation (ELEMENT USES)
%    mixing H & P, see line 298ff below (or search for ***).
%    Fix by: Sven Leyffer, U. Dundee, September 2000
% 
%    The elastodynamic lubrification problem by Kostreva.
% 
%    Source:
%    M.M. Kostreva,
%    "Elasto-hydrodynamic lubrification: a non-linear
%    complementarity problem",
%    International Journal for Numerical Methods in Fluids,
%    4: 377-397, 1984.
% 
%    This problem is problem #5 in More's test set.
% 
%    SIF input: Ph. Toint, June 1990.
% 
%    classification = 'C-CQOR2-MN-V-V'
% 
%    Number of discretized points per unit length
% 
%       Alternative values for the SIF file parameters:
% IE NN                  10             $-PARAMETER n = 151    original value
% IE NN                  50             $-PARAMETER n = 751
% IE NN                  250            $-PARAMETER n = 3751
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'LUBRIFC';

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
            v_('NN') = 5;  %  SIF file default value
        else
            v_('NN') = varargin{1};
        end
        v_('ALPHA') = 1.838;
        v_('LAMBDA') = 1.642;
        v_('XA') = -3.0;
        v_('XF') = 2.0;
        v_('N') = 5*v_('NN');
        v_('0') = 0;
        v_('1') = 1;
        v_('2') = 2;
        v_('3') = 3;
        v_('4') = 4;
        v_('PI') = 3.1415926535;
        v_('2N') = 2*v_('N');
        v_('2N-2') = -2+v_('2N');
        v_('2N-1') = -1+v_('2N');
        v_('2N+2') = 2+v_('2N');
        v_('-XA') = -1.0*v_('XA');
        v_('LEN') = v_('XF')+v_('-XA');
        v_('1/PI') = 1.0/v_('PI');
        v_('1/2PI') = 0.5*v_('1/PI');
        v_('RN') = v_('N');
        v_('1/N') = 1.0/v_('RN');
        v_('DX') = v_('LEN')*v_('1/N');
        v_('1/DX') = 1.0/v_('DX');
        v_('L/DX') = v_('LAMBDA')*v_('1/DX');
        v_('-L/DX') = -1.0*v_('L/DX');
        v_('1/DX2') = v_('1/DX')*v_('1/DX');
        v_('-1/DX2') = -1.0*v_('1/DX2');
        v_('DX/PI') = v_('DX')*v_('1/PI');
        v_('2DX/PI') = 2.0*v_('DX/PI');
        v_('DX/2') = 0.5*v_('DX');
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        [iv,ix_] = s2mpjlib('ii','K',ix_);
        pb.xnames{iv} = 'K';
        for I=v_('0'):v_('2'):v_('2N')
            [iv,ix_] = s2mpjlib('ii',['P',int2str(I)],ix_);
            pb.xnames{iv} = ['P',int2str(I)];
        end
        for J=v_('1'):v_('2'):v_('2N-1')
            [iv,ix_] = s2mpjlib('ii',['H',int2str(J)],ix_);
            pb.xnames{iv} = ['H',int2str(J)];
        end
        for I=v_('2'):v_('2'):v_('2N-2')
            [iv,ix_] = s2mpjlib('ii',['R',int2str(I)],ix_);
            pb.xnames{iv} = ['R',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for I=v_('2'):v_('2'):v_('2N-2')
            [ig,ig_] = s2mpjlib('ii',['R',int2str(round(v_('0')))],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['R',int2str(round(v_('0')))];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['P',int2str(I)]);
            valA(end+1) = v_('2DX/PI');
        end
        [ig,ig_] = s2mpjlib('ii',['R',int2str(round(v_('0')))],ig_);
        gtype{ig}  = '==';
        cnames{ig} = ['R',int2str(round(v_('0')))];
        irA(end+1)  = ig;
        icA(end+1)  = ix_(['P',int2str(round(v_('2N')))]);
        valA(end+1) = v_('DX/PI');
        [ig,ig_] = s2mpjlib('ii','COMPL',ig_);
        gtype{ig} = '<>';
        for I=v_('2'):v_('2'):v_('2N-2')
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            [ig,ig_] = s2mpjlib('ii',['DR',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['DR',int2str(I)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['H',int2str(round(v_('I+1')))]);
            valA(end+1) = v_('L/DX');
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['H',int2str(round(v_('I-1')))]);
            valA(end+1) = v_('-L/DX');
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['R',int2str(I)]);
            valA(end+1) = -1.0;
        end
        for J=v_('1'):v_('2'):v_('2N-1')
            v_('-J') = -1*J;
            [ig,ig_] = s2mpjlib('ii',['DH',int2str(J)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['DH',int2str(J)];
            irA(end+1)  = ig;
            icA(end+1)  = ix_('K');
            valA(end+1) = 1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_(['H',int2str(J)]);
            valA(end+1) = -1.0;
            for I=v_('2'):v_('2N')
                v_(['C',int2str(I)]) = 0.0;
            end
            v_('RI-J') = v_('-J');
            v_('I-JDX') = v_('RI-J')*v_('DX/2');
            v_('ALN') = abs(v_('I-JDX'));
            v_('LN') = log(v_('ALN'));
            v_('T1') = v_('I-JDX')*v_('LN');
            v_('COEFF') = v_('T1')*v_('1/2PI');
            v_(['C',int2str(round(v_('2')))]) = v_(['C',int2str(round(v_('2')))])+...
                 v_('COEFF');
            v_('I-J') = 2+v_('-J');
            v_('RI-J') = v_('I-J');
            v_('I-JDX') = v_('RI-J')*v_('DX/2');
            v_('ALN') = abs(v_('I-JDX'));
            v_('LN') = log(v_('ALN'));
            v_('T1') = v_('I-JDX')*v_('LN');
            v_('COEFF') = v_('T1')*v_('1/PI');
            v_(['C',int2str(round(v_('4')))]) = v_(['C',int2str(round(v_('4')))])+...
                 v_('COEFF');
            for I=v_('4'):v_('2'):v_('2N-2')
                v_('I-2') = -2+I;
                v_('I+2') = 2+I;
                v_('I-J') = I+v_('-J');
                v_('RI-J') = v_('I-J');
                v_('I-JDX') = v_('RI-J')*v_('DX/2');
                v_('ALN') = abs(v_('I-JDX'));
                v_('LN') = log(v_('ALN'));
                v_('T1') = v_('I-JDX')*v_('LN');
                v_('COEFF') = v_('T1')*v_('1/PI');
                v_(['C',int2str(round(v_('I+2')))]) = v_(['C',int2str(round(v_('I+2')))])+...
                     v_('COEFF');
                v_('-COEFF') = -1.0*v_('COEFF');
                v_(['C',int2str(round(v_('I-2')))]) = v_(['C',int2str(round(v_('I-2')))])+...
                     v_('-COEFF');
            end
            v_('I-J') = v_('2N')+v_('-J');
            v_('RI-J') = v_('I-J');
            v_('I-JDX') = v_('RI-J')*v_('DX/2');
            v_('ALN') = abs(v_('I-JDX'));
            v_('LN') = log(v_('ALN'));
            v_('T1') = v_('I-JDX')*v_('LN');
            v_('COEFF') = v_('T1')*v_('1/2PI');
            v_('-COEFF') = -1.0*v_('COEFF');
            v_(['C',int2str(round(v_('2N-2')))]) = v_(['C',int2str(round(v_('2N-2')))])+...
                 v_('-COEFF');
            for I=v_('2'):v_('2'):v_('2N-2')
                [ig,ig_] = s2mpjlib('ii',['DH',int2str(J)],ig_);
                gtype{ig}  = '==';
                cnames{ig} = ['DH',int2str(J)];
                irA(end+1)  = ig;
                icA(end+1)  = ix_(['P',int2str(I)]);
                valA(end+1) = v_(['C',int2str(I)]);
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
        pbm.gconst(ig_(['R',int2str(round(v_('0')))])) = 1.0;
        for J=v_('1'):v_('2'):v_('2N-1')
            v_('RJ') = J;
            v_('JDX') = v_('RJ')*v_('DX/2');
            v_('XJ') = v_('XA')+v_('JDX');
            v_('XJSQ') = v_('XJ')*v_('XJ');
            v_('XJSQ+1') = 1.0+v_('XJSQ');
            v_('RHS') = -1.0*v_('XJSQ+1');
            pbm.gconst(ig_(['DH',int2str(J)])) = v_('RHS');
        end
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('K')) = -Inf;
        pb.xupper(ix_('K'),1) = +Inf;
        for I=v_('2'):v_('2'):v_('2N-2')
            pb.xupper(ix_(['P',int2str(I)])) = 3.0;
            pb.xlower(ix_(['P',int2str(I)]),1) = 0.0;
        end
        for I=v_('1'):v_('2'):v_('2N-1')
            pb.xlower(ix_(['H',int2str(I)])) = -Inf;
            pb.xupper(ix_(['H',int2str(I)]),1) = +Inf;
        end
        pb.xlower(ix_(['P',int2str(round(v_('0')))]),1) = 0.0;
        pb.xupper(ix_(['P',int2str(round(v_('0')))]),1) = 0.0;
        pb.xlower(ix_(['P',int2str(round(v_('2N')))]),1) = 0.0;
        pb.xupper(ix_(['P',int2str(round(v_('2N')))]),1) = 0.0;
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 0.0*ones(pb.n,1);
        v_('2NN') = v_('NN')+v_('NN');
        v_('4NN') = 4*v_('NN');
        for I=v_('2'):v_('2'):v_('4NN')
            v_('RI') = I;
            v_('IDX') = v_('RI')*v_('DX/2');
            v_('XI') = v_('XA')+v_('IDX');
            v_('LIN') = 0.02*v_('XI');
            v_('PI0') = 0.06+v_('LIN');
            if(isKey(ix_,['P',int2str(I)]))
                pb.x0(ix_(['P',int2str(I)]),1) = v_('PI0');
            else
                pb.y0(find(pbm.congrps==ig_(['P',int2str(I)])),1) = v_('PI0');
            end
        end
        v_('4NN+2') = 2+v_('4NN');
        v_('8NN') = 8*v_('NN');
        for I=v_('4NN+2'):v_('2'):v_('8NN')
            v_('RI') = I;
            v_('IDX') = v_('RI')*v_('DX/2');
            v_('XI') = v_('XA')+v_('IDX');
            v_('XISQ') = v_('XI')*v_('XI');
            v_('-XISQ') = -1.0*v_('XISQ');
            v_('1-XISQ') = 1.0+v_('-XISQ');
            v_('PI0') = sqrt(v_('1-XISQ'));
            if(isKey(ix_,['P',int2str(I)]))
                pb.x0(ix_(['P',int2str(I)]),1) = v_('PI0');
            else
                pb.y0(find(pbm.congrps==ig_(['P',int2str(I)])),1) = v_('PI0');
            end
        end
        v_('8NN+2') = 2+v_('8NN');
        for I=v_('8NN+2'):v_('2'):v_('2N')
            if(isKey(ix_,['P',int2str(I)]))
                pb.x0(ix_(['P',int2str(I)]),1) = 0.0;
            else
                pb.y0(find(pbm.congrps==ig_(['P',int2str(I)])),1) = 0.0;
            end
        end
        for J=v_('1'):v_('2'):v_('2N-1')
            v_('RJ') = J;
            v_('JDX') = v_('RJ')*v_('DX/2');
            v_('XJ') = v_('XA')+v_('JDX');
            v_('XJSQ') = v_('XJ')*v_('XJ');
            if(isKey(ix_,['H',int2str(J)]))
                pb.x0(ix_(['H',int2str(J)]),1) = v_('XJSQ');
            else
                pb.y0(find(pbm.congrps==ig_(['H',int2str(J)])),1) = v_('XJSQ');
            end
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eREY',iet_);
        elftv{it}{1} = 'PA';
        elftv{it}{2} = 'PB';
        elftv{it}{3} = 'H';
        elftp{it}{1} = 'A';
        [it,iet_] = s2mpjlib( 'ii', 'en2PR',iet_);
        elftv{it}{1} = 'P';
        elftv{it}{2} = 'R';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for J=v_('1'):v_('2'):v_('2N-1')
            v_('I+') = 1+J;
            v_('I-') = -1+J;
            ename = ['ER',int2str(J)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eREY';
            ielftype(ie) = iet_('eREY');
            vname = ['P',int2str(round(v_('I-')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
            posev = find(strcmp('PA',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['H',int2str(J)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
            posev = find(strcmp('H',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['P',int2str(round(v_('I+')))];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
            posev = find(strcmp('PB',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('A',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('ALPHA');
        end
        for I=v_('2'):v_('2'):v_('2N-2')
            ename = ['EC',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'en2PR';
            ielftype(ie) = iet_('en2PR');
            vname = ['P',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
            posev = find(strcmp('P',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = ['R',int2str(I)];
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],0.0);
            posev = find(strcmp('R',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('2'):v_('2'):v_('2N-2')
            ig = ig_('COMPL');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['EC',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        for I=v_('2'):v_('2'):v_('2N-2')
            v_('I-1') = -1+I;
            v_('I+1') = 1+I;
            ig = ig_(['DR',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['ER',int2str(round(v_('I-1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('1/DX2');
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['ER',int2str(round(v_('I+1')))]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = v_('-1/DX2');
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
        pb.objlower = 0.0;
%    Solution
% LO SOLTN                0.0
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CQOR2-MN-V-V';
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

    case 'eREY'

        EV_  = varargin{1};
        iel_ = varargin{2};
        HA = -0.5*pbm.elpar{iel_}(1);
        EARG = HA*(EV_(1)+EV_(2));
        E = exp(EARG);
        PAMPB = EV_(1)-EV_(2);
        T1 = PAMPB*HA+1.0;
        T2 = PAMPB*HA-1.0;
        HSQ = EV_(3)*EV_(3);
        HCB = HSQ*EV_(3);
        varargout{1} = PAMPB*HCB*E;
        if(nargout>1)
            g_(1,1) = T1*HCB*E;
            g_(2,1) = T2*HCB*E;
            g_(3,1) = 3.0*PAMPB*HSQ*E;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,1) = HCB*E*HA*(T1+1.0);
                H_(1,2) = HCB*E*HA*(T1-1.0);
                H_(2,1) = H_(1,2);
                H_(1,3) = 3.0*T1*HSQ*E;
                H_(3,1) = H_(1,3);
                H_(2,2) = HCB*E*HA*(T2-1.0);
                H_(2,3) = 3.0*T2*HSQ*E;
                H_(3,2) = H_(2,3);
                H_(3,3) = 6.0*EV_(3)*PAMPB*E;
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

