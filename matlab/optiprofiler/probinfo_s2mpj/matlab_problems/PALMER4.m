function varargout = PALMER4(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : PALMER4
%    *********
% 
%    A nonlinear least squares problem with bounds
%    arising from chemical kinetics.
% 
%    model: H-N=C=Se TZVP + MP2
%    fitting Y to A X**2 + B / ( C + X**2 / D ), B, C, D nonnegative.
% 
%    Source:
%    M. Palmer, Edinburgh, private communication.
% 
%    SIF input: Nick Gould, 1990.
% 
%    classification = 'SBR2-RN-4-0'
% 
%    Number of data points
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'PALMER4';

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
        v_  = configureDictionary('string','double');
        ix_ = configureDictionary('string','double');
        ig_ = configureDictionary('string','double');
        v_('M') = 23;
        v_('1') = 1;
        v_('X1') = -1.658063;
        v_('X2') = -1.570796;
        v_('X3') = -1.396263;
        v_('X4') = -1.221730;
        v_('X5') = -1.047198;
        v_('X6') = -0.872665;
        v_('X7') = -0.741119;
        v_('X8') = -0.698132;
        v_('X9') = -0.523599;
        v_('X10') = -0.349066;
        v_('X11') = -0.174533;
        v_('X12') = 0.0;
        v_('X13') = 0.174533;
        v_('X14') = 0.349066;
        v_('X15') = 0.523599;
        v_('X16') = 0.698132;
        v_('X17') = 0.741119;
        v_('X18') = 0.872665;
        v_('X19') = 1.047198;
        v_('X20') = 1.221730;
        v_('X21') = 1.396263;
        v_('X22') = 1.570796;
        v_('X23') = 1.658063;
        v_('Y1') = 67.27625;
        v_('Y2') = 52.8537;
        v_('Y3') = 30.2718;
        v_('Y4') = 14.9888;
        v_('Y5') = 5.5675;
        v_('Y6') = 0.92603;
        v_('Y7') = 0.0;
        v_('Y8') = 0.085108;
        v_('Y9') = 1.867422;
        v_('Y10') = 5.014768;
        v_('Y11') = 8.263520;
        v_('Y12') = 9.8046208;
        v_('Y13') = 8.263520;
        v_('Y14') = 5.014768;
        v_('Y15') = 1.867422;
        v_('Y16') = 0.085108;
        v_('Y17') = 0.0;
        v_('Y18') = 0.92603;
        v_('Y19') = 5.5675;
        v_('Y20') = 14.9888;
        v_('Y21') = 30.2718;
        v_('Y22') = 52.8537;
        v_('Y23') = 67.27625;
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
        pb.n   = numEntries(ix_);
        ngrp   = numEntries(ig_);
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
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eQUOT',iet_);
        elftv{it}{1} = 'B';
        elftv{it}{2} = 'C';
        elftv{it}{3} = 'D';
        elftp{it}{1} = 'XSQR';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
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
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
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
% LO SOLTN               2285.383227
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'SBR2-RN-4-0';
        %%%%%%%%%%% REDUCED-PRECISION CONVERSION %%%%%%%%%%%
        if(strcmp(action,'setup_redprec'))
            varargout{1} = s2mpjlib('convert',pb,  pbm.ndigs);
            varargout{2} = s2mpjlib('convert',pbm, pbm.ndigs);
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
          'cJxv','Lxy','Lgxy','LgHxy','LIxy','LIgxy','LIgHxy','LHxyv','LIHxyv'}

        if(isfield(pbm,'name')&&strcmp(pbm.name,name))
            pbm.has_globs = [0,0];
            [varargout{1:max(1,nargout)}] = s2mpjlib(action,pbm,varargin{:});
        else
            disp(['ERROR: please run ',name,' with action = setup'])
            [varargout{1:nargout}] = deal(repmat(NaN,1:nargout));
        end

    otherwise
        disp([' ERROR: unknown action ',action,' requested from ',name,'.m'])
    end

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

