function varargout = PALMER2E(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : PALMER2E
%    *********
% 
%    A nonlinear least squares problem
%    arising from chemical kinetics.
% 
%     model: H-N=C=O TZVP + MP2
%    fitting Y to A0 + A2 X**2 + A4 X**4 + A6 X**6 + A8 X**8 +
%                 A10 X**10 + L * EXP( -K X**2 )
% 
%    Source:
%    M. Palmer, Edinburgh, private communication.
% 
%    SIF input: Nick Gould, 1990.
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

name = 'PALMER2E';

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
        v_('X1') = -1.745329;
        v_('X2') = -1.570796;
        v_('X3') = -1.396263;
        v_('X4') = -1.221730;
        v_('X5') = -1.047198;
        v_('X6') = -0.937187;
        v_('X7') = -0.872665;
        v_('X8') = -0.698132;
        v_('X9') = -0.523599;
        v_('X10') = -0.349066;
        v_('X11') = -0.174533;
        v_('X12') = 0.0;
        v_('X13') = 0.174533;
        v_('X14') = 0.349066;
        v_('X15') = 0.523599;
        v_('X16') = 0.698132;
        v_('X17') = 0.872665;
        v_('X18') = 0.937187;
        v_('X19') = 1.047198;
        v_('X20') = 1.221730;
        v_('X21') = 1.396263;
        v_('X22') = 1.570796;
        v_('X23') = 1.745329;
        v_('Y1') = 72.676767;
        v_('Y2') = 40.149455;
        v_('Y3') = 18.8548;
        v_('Y4') = 6.4762;
        v_('Y5') = 0.8596;
        v_('Y6') = 0.00000;
        v_('Y7') = 0.2730;
        v_('Y8') = 3.2043;
        v_('Y9') = 8.1080;
        v_('Y10') = 13.4291;
        v_('Y11') = 17.7149;
        v_('Y12') = 19.4529;
        v_('Y13') = 17.7149;
        v_('Y14') = 13.4291;
        v_('Y15') = 8.1080;
        v_('Y16') = 3.2053;
        v_('Y17') = 0.2730;
        v_('Y18') = 0.00000;
        v_('Y19') = 0.8596;
        v_('Y20') = 6.4762;
        v_('Y21') = 18.8548;
        v_('Y22') = 40.149455;
        v_('Y23') = 72.676767;
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
        for I=v_('1'):v_('M')
            v_('XSQR') = v_(['X',int2str(I)])*v_(['X',int2str(I)]);
            v_('XQUART') = v_('XSQR')*v_('XSQR');
            v_('X**6') = v_('XSQR')*v_('XQUART');
            v_('X**8') = v_('XSQR')*v_('X**6');
            v_('X**10') = v_('XSQR')*v_('X**8');
            v_('X**12') = v_('XSQR')*v_('X**10');
            v_('X**14') = v_('XSQR')*v_('X**12');
            [ig,ig_] = s2mpjlib('ii',['O',int2str(I)],ig_);
            gtype{ig} = '<>';
            irA(end+1)  = ig;
            icA(end+1)  = ix_('A0');
            valA(end+1) = 1.0;
            irA(end+1)  = ig;
            icA(end+1)  = ix_('A2');
            valA(end+1) = v_('XSQR');
            irA(end+1)  = ig;
            icA(end+1)  = ix_('A4');
            valA(end+1) = v_('XQUART');
            irA(end+1)  = ig;
            icA(end+1)  = ix_('A6');
            valA(end+1) = v_('X**6');
            irA(end+1)  = ig;
            icA(end+1)  = ix_('A8');
            valA(end+1) = v_('X**8');
            irA(end+1)  = ig;
            icA(end+1)  = ix_('A10');
            valA(end+1) = v_('X**10');
        end
        %%%%%%%%%%%%%%% GLOBAL DIMENSIONS %%%%%%%%%%%%%%%%%
        pb.n   = ix_.Count;
        ngrp   = ig_.Count;
        pbm.objgrps = [1:ngrp];
        pb.m        = 0;
        %%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%
        pbm.gconst = zeros(ngrp,1);
        for I=v_('1'):v_('M')
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
        for I=v_('1'):v_('M')
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
        for I=v_('1'):v_('M')
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
% LO SOLTN                2.0650351D-04
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

