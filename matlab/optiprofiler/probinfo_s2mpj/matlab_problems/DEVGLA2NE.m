function varargout = DEVGLA2NE(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : DEVGLA2NE
%    *********
% 
%    SCIPY global optimization benchmark example DeVilliersGlasser02
% 
%    Fit: y  = x_1 x_2^t tanh ( t x_3 + sin( t x_4 ) ) cos( t e^x_5 )  +  e
% 
%    Source:  Problem from the SCIPY benchmark set
%      https://github.com/scipy/scipy/tree/master/benchmarks/ ...
%              benchmarks/go_benchmark_functions
% 
%    Nonlinear-equation formulation of DEVGLA2.SIF
% 
%    SIF input: Nick Gould, Jan 2020
% 
%    classification = 'C-CNOR2-MN-5-16'
% 
%    Number of data values
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'DEVGLA2NE';

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
        v_('M') = 16;
        v_('N') = 5;
        v_('1') = 1;
        v_('A') = 1.27;
        v_('LNA') = log(v_('A'));
        for I=v_('1'):v_('M')
            v_('RI') = I;
            v_('RIM1') = -1.0+v_('RI');
            v_('T') = 0.1*v_('RIM1');
            v_(['T',int2str(I)]) = v_('T');
            v_('TLNA') = v_('T')*v_('LNA');
            v_('AT') = exp(v_('TLNA'));
            v_('TP') = 3.012*v_('T');
            v_('TP2') = 2.13*v_('T');
            v_('STP2') = sin(v_('TP2'));
            v_('TPA') = v_('TP')+v_('STP2');
            v_('HTPA') = tanh(v_('TPA'));
            v_('EC') = exp(0.507);
            v_('ECT') = v_('EC')*v_('T');
            v_('CECT') = cos(v_('ECT'));
            v_('P') = v_('AT')*v_('HTPA');
            v_('PP') = v_('P')*v_('CECT');
            v_('PPP') = 53.81*v_('PP');
            v_(['Y',int2str(I)]) = v_('PPP');
        end
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
        for I=v_('1'):v_('M')
            [ig,ig_] = s2mpjlib('ii',['F',int2str(I)],ig_);
            gtype{ig}  = '==';
            cnames{ig} = ['F',int2str(I)];
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
        for I=v_('1'):v_('M')
            pbm.gconst(ig_(['F',int2str(I)])) = v_(['Y',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.y0 = zeros(pb.m,1);
        if(isKey(ix_,'X1'))
            pb.x0(ix_('X1'),1) = 20.0;
        else
            pb.y0(find(pbm.congrps==ig_('X1')),1) = 20.0;
        end
        if(isKey(ix_,'X2'))
            pb.x0(ix_('X2'),1) = 2.0;
        else
            pb.y0(find(pbm.congrps==ig_('X2')),1) = 2.0;
        end
        if(isKey(ix_,'X3'))
            pb.x0(ix_('X3'),1) = 2.0;
        else
            pb.y0(find(pbm.congrps==ig_('X3')),1) = 2.0;
        end
        if(isKey(ix_,'X4'))
            pb.x0(ix_('X4'),1) = 2.0;
        else
            pb.y0(find(pbm.congrps==ig_('X4')),1) = 2.0;
        end
        if(isKey(ix_,'X5'))
            pb.x0(ix_('X5'),1) = 0.2;
        else
            pb.y0(find(pbm.congrps==ig_('X5')),1) = 0.2;
        end
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eDG2',iet_);
        elftv{it}{1} = 'X1';
        elftv{it}{2} = 'X2';
        elftv{it}{3} = 'X3';
        elftv{it}{4} = 'X4';
        elftv{it}{5} = 'X5';
        elftp{it}{1} = 'T';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('M')
            ename = ['E',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eDG2';
            ielftype(ie) = iet_('eDG2');
            vname = 'X1';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X2';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X3',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X4';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X4',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'X5';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('X5',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('T',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['T',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('1'):v_('M')
            ig = ig_(['F',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(I)]);
            nlc = union(nlc,ig);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Least square problems are bounded below by zero
        pb.objlower = 0.0;
%    Solution
% LO SOLUTION            0.0
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%%%%%%%%%% FORM clower AND cupper %%%%%%%%%%%%%
        pb.clower(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        pb.cupper(pb.nle+1:pb.nle+pb.neq) = zeros(pb.neq,1);
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        [~,pb.lincons]  = ismember(setdiff(pbm.congrps,nlc),pbm.congrps);
        pb.pbclass = 'C-CNOR2-MN-5-16';
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

    case 'eDG2'

        EV_  = varargin{1};
        iel_ = varargin{2};
        X2T = EV_(2)^pbm.elpar{iel_}(1);
        F2 = X2T;
        F2X2 = pbm.elpar{iel_}(1)*EV_(2)^(pbm.elpar{iel_}(1)-1.0e0);
        F2X2X2 =...
              pbm.elpar{iel_}(1)*(pbm.elpar{iel_}(1)-1.0e0)*EV_(2)^(pbm.elpar{iel_}(1)-2.0e0);
        X3T = EV_(3)*pbm.elpar{iel_}(1);
        X4T = EV_(4)*pbm.elpar{iel_}(1);
        SINX4T = sin(X4T);
        COSX4T = cos(X4T);
        A = X3T+SINX4T;
        AX3 = pbm.elpar{iel_}(1);
        AX4 = pbm.elpar{iel_}(1)*COSX4T;
        AX4X4 = -pbm.elpar{iel_}(1)*pbm.elpar{iel_}(1)*SINX4T;
        F = tanh(A);
        FA = 1.0/(cosh(A))^2;
        FAA = -2.0*FA*F;
        F3 = F;
        F3X3 = FA*AX3;
        F3X4 = FA*AX4;
        F3X3X3 = FAA*AX3*AX3;
        F3X3X4 = FAA*AX3*AX4;
        F3X4X4 = FA*AX4X4+FAA*AX4*AX4;
        EX5 = exp(EV_(5));
        TEX5 = pbm.elpar{iel_}(1)*EX5;
        STEX5 = sin(TEX5);
        CTEX5 = cos(TEX5);
        F4 = CTEX5;
        F4X5 = -STEX5*TEX5;
        F4X5X5 = -STEX5*TEX5-CTEX5*TEX5*TEX5;
        varargout{1} = EV_(1)*F2*F3*F4;
        if(nargout>1)
            g_(1,1) = F2*F3*F4;
            g_(2,1) = EV_(1)*F2X2*F3*F4;
            g_(3,1) = EV_(1)*F2*F3X3*F4;
            g_(4,1) = EV_(1)*F2*F3X4*F4;
            g_(5,1) = EV_(1)*F2*F3*F4X5;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(5,5);
                H_(1,2) = F2X2*F3*F4;
                H_(2,1) = H_(1,2);
                H_(1,3) = F2*F3X3*F4;
                H_(3,1) = H_(1,3);
                H_(1,4) = F2*F3X4*F4;
                H_(4,1) = H_(1,4);
                H_(1,5) = F2*F3*F4X5;
                H_(5,1) = H_(1,5);
                H_(2,2) = EV_(1)*F2X2X2*F3*F4;
                H_(2,3) = EV_(1)*F2X2*F3X3*F4;
                H_(3,2) = H_(2,3);
                H_(2,4) = EV_(1)*F2X2*F3X4*F4;
                H_(4,2) = H_(2,4);
                H_(2,5) = EV_(1)*F2X2*F3*F4X5;
                H_(5,2) = H_(2,5);
                H_(3,3) = EV_(1)*F2*F3X3X3*F4;
                H_(3,4) = EV_(1)*F2*F3X3X4*F4;
                H_(4,3) = H_(3,4);
                H_(3,5) = EV_(1)*F2*F3X3*F4X5;
                H_(5,3) = H_(3,5);
                H_(4,4) = EV_(1)*F2*F3X4X4*F4;
                H_(4,5) = EV_(1)*F2*F3X4*F4X5;
                H_(5,4) = H_(4,5);
                H_(5,5) = EV_(1)*F2*F3*F4X5X5;
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

