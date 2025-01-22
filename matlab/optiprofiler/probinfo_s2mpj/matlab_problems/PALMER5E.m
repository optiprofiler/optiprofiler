function varargout = PALMER5E(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : PALMER5E
%    *********
% 
%    A nonlinear least squares problem
%    arising from chemical kinetics.
% 
%    model: H-N=C=Se TZVP + MP2
%    fitting Y to A0 T_0 + A2 T_2 + A4 T_4 + A6 T_6 + A8 T_8 +
%                 A10 T_10 + A12 T_12 + A14 T_14
%                 + L * EXP( -K X**2 )
%    where T_i is the i-th (shifted) Chebyshev polynomial
% 
%    Source:
%    M.  Palmer, Edinburgh, private communication.
% 
%    SIF input: Nick Gould, 1992.
% 
%    classification = 'C-CSBR2-RN-8-0'
% 
%    Number of data points
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 25 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'PALMER5E';

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
        v_('M') = 23;
        v_('1') = 1;
        v_('2') = 2;
        v_('12') = 12;
        v_('14') = 14;
        v_('X12') = 0.000000;
        v_('X13') = 1.570796;
        v_('X14') = 1.396263;
        v_('X15') = 1.308997;
        v_('X16') = 1.221730;
        v_('X17') = 1.125835;
        v_('X18') = 1.047198;
        v_('X19') = 0.872665;
        v_('X20') = 0.698132;
        v_('X21') = 0.523599;
        v_('X22') = 0.349066;
        v_('X23') = 0.174533;
        v_('B') = v_('X13');
        v_('A') = -1.0e+0*v_('B');
        v_('DIFF') = 2.0e+0*v_('B');
        v_('Y12') = 83.57418;
        v_('Y13') = 81.007654;
        v_('Y14') = 18.983286;
        v_('Y15') = 8.051067;
        v_('Y16') = 2.044762;
        v_('Y17') = 0.000000;
        v_('Y18') = 1.170451;
        v_('Y19') = 10.479881;
        v_('Y20') = 25.785001;
        v_('Y21') = 44.126844;
        v_('Y22') = 62.822177;
        v_('Y23') = 77.719674;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        irA  = [];
        icA  = [];
        valA = [];
        [iv,ix_] = s2mpjlib('ii','A0',ix_);
        pb.xnames{iv} = 'A0';
        [iv,ix_] = s2mpjlib('ii','A2',ix_);
        pb.xnames{iv} = 'A2';
        [iv,ix_] = s2mpjlib('ii','A4',ix_);
        pb.xnames{iv} = 'A4';
        [iv,ix_] = s2mpjlib('ii','A6',ix_);
        pb.xnames{iv} = 'A6';
        [iv,ix_] = s2mpjlib('ii','A8',ix_);
        pb.xnames{iv} = 'A8';
        [iv,ix_] = s2mpjlib('ii','A10',ix_);
        pb.xnames{iv} = 'A10';
        [iv,ix_] = s2mpjlib('ii','K',ix_);
        pb.xnames{iv} = 'K';
        [iv,ix_] = s2mpjlib('ii','L',ix_);
        pb.xnames{iv} = 'L';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        for I=v_('12'):v_('M')
            v_('T0') = 1.0e+0;
            v_('Y') = 2.0e+0*v_(['X',int2str(I)]);
            v_('Y') = v_('Y')-v_('A');
            v_('Y') = v_('Y')-v_('B');
            v_('Y') = v_('Y')/v_('DIFF');
            v_('T1') = v_('Y');
            v_('2Y') = 2.0e+0*v_('Y');
            for J=v_('2'):v_('14')
                v_('J-1') = -1+J;
                v_('J-2') = -2+J;
                v_(['T',int2str(J)]) = v_('2Y')*v_(['T',int2str(round(v_('J-1')))]);
                v_(['T',int2str(J)]) =...
                      v_(['T',int2str(J)])-v_(['T',int2str(round(v_('J-2')))]);
            end
            [ig,ig_] = s2mpjlib('ii',['O',int2str(I)],ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_('A0');
            valA(end+1) = v_('T0');
            irA(end+1)  = ig;
            icA(end+1)  = ix_('A2');
            valA(end+1) = v_('T2');
            irA(end+1)  = ig;
            icA(end+1)  = ix_('A4');
            valA(end+1) = v_('T4');
            irA(end+1)  = ig;
            icA(end+1)  = ix_('A6');
            valA(end+1) = v_('T6');
            irA(end+1)  = ig;
            icA(end+1)  = ix_('A8');
            valA(end+1) = v_('T8');
            irA(end+1)  = ig;
            icA(end+1)  = ix_('A10');
            valA(end+1) = v_('T10');
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('12'):v_('M')
            pbm.gconst(ig_(['O',int2str(I)])) = v_(['Y',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = zeros(pb.n,1);
        pb.xupper = Inf*ones(pb.n,1);
        pb.xlower(ix_('A0')) = -Inf;
        pb.xupper(ix_('A0'),1) = +Inf;
        pb.xlower(ix_('A2')) = -Inf;
        pb.xupper(ix_('A2'),1) = +Inf;
        pb.xlower(ix_('A4')) = -Inf;
        pb.xupper(ix_('A4'),1) = +Inf;
        pb.xlower(ix_('A6')) = -Inf;
        pb.xupper(ix_('A6'),1) = +Inf;
        pb.xlower(ix_('A8')) = -Inf;
        pb.xupper(ix_('A8'),1) = +Inf;
        pb.xlower(ix_('A10')) = -Inf;
        pb.xupper(ix_('A10'),1) = +Inf;
        pb.xlower(ix_('L')) = -Inf;
        pb.xupper(ix_('L'),1) = +Inf;
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 1.0*ones(pb.n,1);
        pb.x0(ix_('A0'),1) = 1.9264e+01;
        pb.x0(ix_('A2'),1) = -1.7302;
        pb.x0(ix_('A4'),1) = 4.0794e+01;
        pb.x0(ix_('A6'),1) = 8.3021e-01;
        pb.x0(ix_('A8'),1) = 3.7090;
        pb.x0(ix_('A10'),1) = -1.7723e-01;
        pb.x0(ix_('K'),1) = 10.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'ePROD',iet_);
        elftv{it}{1} = 'K';
        elftv{it}{2} = 'L';
        elftp{it}{1} = 'XSQR';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = containers.Map('KeyType','char','ValueType','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('12'):v_('M')
            v_('XSQR') = v_(['X',int2str(I)])*v_(['X',int2str(I)]);
            ename = ['E',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'ePROD';
            ielftype(ie) = iet_('ePROD');
            vname = 'K';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('K',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'L';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('L',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('XSQR',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_('XSQR');
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = containers.Map('KeyType','char','ValueType','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal([]);
        nlc = [];
        for I=v_('12'):v_('M')
            ig = ig_(['O',int2str(I)]);
            pbm.grftype{ig} = 'gL2';
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Least square problems are bounded below by zero
        pb.objlower = 0.0;
%    Solution
% LO SOLTN              1.48003482D-04
        %%%%%%%%% BUILD THE SPARSE MATRICES %%%%%%%%%%%%%%%
        pbm.A = sparse(irA,icA,valA,ngrp,pb.n);
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'C-CSBR2-RN-8-0';
        pbm.objderlvl = 2;
        pb.objderlvl = pbm.objderlvl;
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

    case 'ePROD'

        EV_  = varargin{1};
        iel_ = varargin{2};
        EXPON = exp(-EV_(1)*pbm.elpar{iel_}(1));
        varargout{1} = EV_(2)*EXPON;
        if(nargout>1)
            g_(1,1) = -pbm.elpar{iel_}(1)*EV_(2)*EXPON;
            g_(2,1) = EXPON;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = pbm.elpar{iel_}(1)*pbm.elpar{iel_}(1)*EV_(2)*EXPON;
                H_(1,2) = -pbm.elpar{iel_}(1)*EXPON;
                H_(2,1) = H_(1,2);
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

