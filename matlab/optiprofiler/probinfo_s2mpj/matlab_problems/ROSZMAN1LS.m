function varargout = ROSZMAN1LS(action,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%    Problem : ROSZMAN1LS
%    *********
% 
%    NIST Data fitting problem ROSZMAN1.
% 
%    Fit: y =  b1 - b2*x - arctan[b3/(x-b4)]/pi + e
% 
%    Source:  Problem from the NIST nonlinear regression test set
%      http://www.itl.nist.gov/div898/strd/nls/nls_main.shtml
% 
%    SIF input: Nick Gould and Tyrone Rees, Oct 2015
% 
%   Reference: Roszman, L., NIST (1979).  
%     Quantum Defects for Sulfur I Atom.
% 
%    classification = 'SUR2-MN-4-0'
% 
%    Number of data values
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent pbm;

name = 'ROSZMAN1LS';

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
        v_('M') = 25;
        v_('N') = 4;
        v_('1') = 1;
        v_('X1') = -4868.68;
        v_('X2') = -4868.09;
        v_('X3') = -4867.41;
        v_('X4') = -3375.19;
        v_('X5') = -3373.14;
        v_('X6') = -3372.03;
        v_('X7') = -2473.74;
        v_('X8') = -2472.35;
        v_('X9') = -2469.45;
        v_('X10') = -1894.65;
        v_('X11') = -1893.40;
        v_('X12') = -1497.24;
        v_('X13') = -1495.85;
        v_('X14') = -1493.41;
        v_('X15') = -1208.68;
        v_('X16') = -1206.18;
        v_('X17') = -1206.04;
        v_('X18') = -997.92;
        v_('X19') = -996.61;
        v_('X20') = -996.31;
        v_('X21') = -834.94;
        v_('X22') = -834.66;
        v_('X23') = -710.03;
        v_('X24') = -530.16;
        v_('X25') = -464.17;
        v_('Y1') = 0.252429;
        v_('Y2') = 0.252141;
        v_('Y3') = 0.251809;
        v_('Y4') = 0.297989;
        v_('Y5') = 0.296257;
        v_('Y6') = 0.295319;
        v_('Y7') = 0.339603;
        v_('Y8') = 0.337731;
        v_('Y9') = 0.333820;
        v_('Y10') = 0.389510;
        v_('Y11') = 0.386998;
        v_('Y12') = 0.438864;
        v_('Y13') = 0.434887;
        v_('Y14') = 0.427893;
        v_('Y15') = 0.471568;
        v_('Y16') = 0.461699;
        v_('Y17') = 0.461144;
        v_('Y18') = 0.513532;
        v_('Y19') = 0.506641;
        v_('Y20') = 0.505062;
        v_('Y21') = 0.535648;
        v_('Y22') = 0.533726;
        v_('Y23') = 0.568064;
        v_('Y24') = 0.612886;
        v_('Y25') = 0.624169;
        %%%%%%%%%%%%%%%%%%%%  VARIABLES %%%%%%%%%%%%%%%%%%%%
        pb.xnames = {};
        for I=v_('1'):v_('N')
            [iv,ix_] = s2mpjlib('ii',['B',int2str(I)],ix_);
            pb.xnames{iv} = ['B',int2str(I)];
        end
        %%%%%%%%%%%%%%%%%%%  DATA GROUPS %%%%%%%%%%%%%%%%%%%
        pbm.A = sparse(0,0);
        for I=v_('1'):v_('M')
            [ig,ig_] = s2mpjlib('ii',['F',int2str(I)],ig_);
            gtype{ig} = '<>';
            iv = ix_('B1');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = 1.0+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = 1.0;
            end
            v_('-X') = -1.0*v_(['X',int2str(I)]);
            iv = ix_('B2');
            if(size(pbm.A,1)>=ig&&size(pbm.A,2)>=iv)
                pbm.A(ig,iv) = v_('-X')+pbm.A(ig,iv);
            else
                pbm.A(ig,iv) = v_('-X');
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
            pbm.gconst(ig_(['F',int2str(I)])) = v_(['Y',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%  BOUNDS %%%%%%%%%%%%%%%%%%%%%
        pb.xlower = -Inf*ones(pb.n,1);
        pb.xupper = +Inf*ones(pb.n,1);
        pb.xlower = zeros(pb.n,1);
        %%%%%%%%%%%%%%%%%%%% START POINT %%%%%%%%%%%%%%%%%%
        pb.x0(1:pb.n,1) = zeros(pb.n,1);
        pb.x0(ix_('B1'),1) = 0.1;
        pb.x0(ix_('B2'),1) = -0.00001;
        pb.x0(ix_('B3'),1) = 1000.0;
        pb.x0(ix_('B4'),1) = -100.0;
        %%%%%%%%%%%%%%%%%%%%% ELFTYPE %%%%%%%%%%%%%%%%%%%%%
        iet_ = configureDictionary('string','double');
        [it,iet_] = s2mpjlib( 'ii', 'eE7',iet_);
        elftv{it}{1} = 'V1';
        elftv{it}{2} = 'V2';
        elftp{it}{1} = 'X';
        %%%%%%%%%%%%%%%%%%% ELEMENT USES %%%%%%%%%%%%%%%%%%
        ie_ = configureDictionary('string','double');
        pbm.elftype = {};
        ielftype    = [];
        pbm.elvar   = {};
        pbm.elpar   = {};
        for I=v_('1'):v_('M')
            ename = ['E',int2str(I)];
            [ie,ie_] = s2mpjlib('ii',ename,ie_);
            pbm.elftype{ie} = 'eE7';
            ielftype(ie) = iet_('eE7');
            vname = 'B3';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V1',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            vname = 'B4';
            [iv,ix_,pb] = s2mpjlib('nlx',vname,ix_,pb,1,[],[],[]);
            posev = find(strcmp('V2',elftv{ielftype(ie)}));
            pbm.elvar{ie}(posev) = iv;
            [~,posep] = ismember('X',elftp{ielftype(ie)});
            pbm.elpar{ie}(posep) = v_(['X',int2str(I)]);
        end
        %%%%%%%%%%%%%%%%%%%%%% GRFTYPE %%%%%%%%%%%%%%%%%%%%
        igt_ = configureDictionary('string','double');
        [it,igt_] = s2mpjlib('ii','gL2',igt_);
        %%%%%%%%%%%%%%%%%%%% GROUP USES %%%%%%%%%%%%%%%%%%%
        [pbm.grelt{1:ngrp}] = deal(repmat([],1,ngrp));
        nlc = [];
        for ig = 1:ngrp
            pbm.grftype{ig} = 'gL2';
        end
        for I=v_('1'):v_('M')
            ig = ig_(['F',int2str(I)]);
            posel = length(pbm.grelt{ig})+1;
            pbm.grelt{ig}(posel) = ie_(['E',int2str(I)]);
            pbm.grelw{ig}(posel) = 1.;
        end
        %%%%%%%%%%%%%%%%%%% OBJECT BOUNDS %%%%%%%%%%%%%%%%%
%    Least square problems are bounded below by zero
        pb.objlower = 0.0;
%    Solution
% LO SOLTN               
        %%%%%%%%% DEFAULT FOR MISSING SECTION(S) %%%%%%%%%%
        %%%%%% RETURN VALUES FROM THE SETUP ACTION %%%%%%%%
        pb.pbclass = 'SUR2-MN-4-0';
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

    case 'e_globs'

        pbm = varargin{1};
        pbm.efpar = [];
        pbm.efpar(1) = 4.0*atan(1.0e0);
        varargout{1} = pbm;

    case 'eE7'

        EV_  = varargin{1};
        iel_ = varargin{2};
        V12 = EV_(1)*EV_(1);
        V13 = EV_(1)*V12;
        V2MX = EV_(2)-pbm.elpar{iel_}(1);
        V2MX2 = V2MX*V2MX;
        V2MX3 = V2MX*V2MX2;
        R = V12/V2MX2+1.0;
        PIR = pbm.efpar(1)*R;
        PIR2 = PIR*R;
        varargout{1} = -atan(EV_(1)/V2MX)/pbm.efpar(1);
        if(nargout>1)
            g_(1,1) = -1.0/(PIR*V2MX);
            g_(2,1) = EV_(1)/(PIR*V2MX2);
            varargout{2} = g_;
            if(nargout>2)
                H_ = sparse(2,2);
                H_(1,1) = 2.0*EV_(1)/(PIR2*V2MX3);
                H_(1,2) = 1.0/(PIR*V2MX2)-2.0*V12/(PIR2*V2MX^4);
                H_(2,1) = H_(1,2);
                H_(2,2) = 2.0*V13/(PIR2*V2MX^5)-2.0*EV_(1)/(PIR*V2MX3);
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
            pbm.has_globs = [1,0];
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

