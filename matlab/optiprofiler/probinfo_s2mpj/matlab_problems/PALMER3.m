function varargout = PALMER3(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : PALMER3
%    *********
% 
%    A nonlinear least squares problem with bounds
%    arising from chemical kinetics.
% 
%     model: H-N=C=S TZVP + MP2
%    fitting Y to A X**2 + B / ( C + X**2 / D ), B, C, D nonnegative.
% 
%    Source:
%    M. Palmer, Edinburgh, private communication.
% 
%    SIF input: Nick Gould, 1990.
% 
%    classification = 'C-CSBR2-RN-4-0'
% 
%    Number of data points
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Translated to Matlab by S2MPJ version 9 XI 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'PALMER3';

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
        v_('M') = 23;
        v_('1') = 1;
        v_('X1') = -1.658063;
        v_('X2') = -1.570796;
        v_('X3') = -1.396263;
        v_('X4') = -1.221730;
        v_('X5') = -1.047198;
        v_('X6') = -0.872665;
        v_('X7') = -0.766531;
        v_('X8') = -0.698132;
        v_('X9') = -0.523599;
        v_('X10') = -0.349066;
        v_('X11') = -0.174533;
        v_('X12') = 0.0;
        v_('X13') = 0.174533;
        v_('X14') = 0.349066;
        v_('X15') = 0.523599;
        v_('X16') = 0.698132;
        v_('X17') = 0.766531;
        v_('X18') = 0.872665;
        v_('X19') = 1.047198;
        v_('X20') = 1.221730;
        v_('X21') = 1.396263;
        v_('X22') = 1.570796;
        v_('X23') = 1.658063;
        v_('Y1') = 64.87939;
        v_('Y2') = 50.46046;
        v_('Y3') = 28.2034;
        v_('Y4') = 13.4575;
        v_('Y5') = 4.6547;
        v_('Y6') = 0.59447;
        v_('Y7') = 0.0000;
        v_('Y8') = 0.2177;
        v_('Y9') = 2.3029;
        v_('Y10') = 5.5191;
        v_('Y11') = 8.5519;
        v_('Y12') = 9.8919;
        v_('Y13') = 8.5519;
        v_('Y14') = 5.5191;
        v_('Y15') = 2.3029;
        v_('Y16') = 0.2177;
        v_('Y17') = 0.0000;
        v_('Y18') = 0.59447;
        v_('Y19') = 4.6547;
        v_('Y20') = 13.4575;
        v_('Y21') = 28.2034;
        v_('Y22') = 50.46046;
        v_('Y23') = 64.87939;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        [iv,ix_] = s2mpjlib('ii','A',ix_);
        pb.xnames{iv} = 'A';
        [iv,ix_] = s2mpjlib('ii','B',ix_);
        pb.xnames{iv} = 'B';
        [iv,ix_] = s2mpjlib('ii','C',ix_);
        pb.xnames{iv} = 'C';
        [iv,ix_] = s2mpjlib('ii','D',ix_);
        pb.xnames{iv} = 'D';
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('M')
            v_('XSQR') = v_(['X',int2str(I)])*v_(['X',int2str(I)]);
            [ig,ig_] = s2mpjlib('ii',['O',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_('A');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('XSQR')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('XSQR');
            end
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
        pb.xlower(ix_('A')) = -Inf;
        pb.xupper(ix_('A'),1) = +Inf;
        pb.xlower(ix_('B'),1) = 0.00001;
        pb.xlower(ix_('C'),1) = 0.00001;
        pb.xlower(ix_('D'),1) = 0.00001;
        %%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0 = 1.0*ones(pb.n,1);
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = containers.Map('KeyType', 'char', 'ValueType','double');
        [it,iet_] = s2mpjlib( 'ii', 'eQUOT',iet_);
        elftv{it}{1} = 'B';
        elftv{it}{2} = 'C';
        elftv{it}{3} = 'D';
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
            pbm.elftype{ie} = 'eQUOT';
            ielftype(ie) = iet_('eQUOT');
            vname = 'B';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('B',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'C';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('C',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'D';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],1.0);
            posev = find(strcmp('D',elftv{ielftype(ie)}));
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
%    SOLTN               2265.95822
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'C-CSBR2-RN-4-0';
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

    case 'eQUOT'

        EV_  = varargin{1};
        iel_ = varargin{2};
        DENOM = 1.0/(EV_(2)+pbm.elpar{iel_}(1)/EV_(3));
        varargout{1} = EV_(1)*DENOM;
        if(nargout>1)
            g_(1,1) = DENOM;
            g_(2,1) = -EV_(1)*DENOM*DENOM;
            g_(3,1) = EV_(1)*pbm.elpar{iel_}(1)*(DENOM/EV_(3))^2;
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(3,3);
                H_(1,2) = -DENOM*DENOM;
                H_(2,1) = H_(1,2);
                H_(1,3) = pbm.elpar{iel_}(1)*(DENOM/EV_(3))^2;
                H_(3,1) = H_(1,3);
                H_(2,2) = 2.0*EV_(1)*DENOM^3;
                H_(2,3) = -2.0*EV_(1)*pbm.elpar{iel_}(1)*DENOM*(DENOM/EV_(3))^2;
                H_(3,2) = H_(2,3);
                H_(3,3) = 2.0*EV_(1)*DENOM^3*(pbm.elpar{iel_}(1)/EV_(3)^2)^2-2.0*...
                     EV_(1)*DENOM^2*pbm.elpar{iel_}(1)/EV_(3)^3;
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

